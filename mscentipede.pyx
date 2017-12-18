import numpy as np
cimport numpy as np
import cython
cimport cython
import cvxopt as cvx
from cvxopt import solvers
from scipy.special import digamma, gammaln, polygamma
import scipy.optimize as spopt
import sys, time, math, pdb

# suppress optimizer output
solvers.options['show_progress'] = False
solvers.options['maxiters'] = 20

# defining some constants
EPS = np.finfo(np.double).tiny
MAX = np.finfo(np.double).max

# defining some simple functions
logistic = lambda x: 1./(1+np.exp(x))
insum = lambda x,axes: np.apply_over_axes(np.sum,x,axes)

def nplog(x):
    """Compute the natural logarithm, handling very
    small floats appropriately.

    """
    try:
        x[x<EPS] = EPS
    except TypeError:
        x = max([x,EPS])
    return np.log(x)


cdef class Data:
    """
    A data structure to store a multiscale representation of
    chromatin accessibility read counts across `N` genomic windows of
    length `L` in `R` replicates.

    Arguments
        reads : array

    """

    def __cinit__(self):

        self.N = 0
        self.L = 0
        self.R = 0
        self.J = 0
        self.valueA = dict()
        self.valueB = dict()
        self.total = dict()

    cdef transform_to_multiscale(self, np.ndarray[np.float64_t, ndim=3] reads):
        """Transform a vector of read counts
        into a multiscale representation.
        
        .. note::
            See msCentipede manual for more details.   

        """

        cdef long k, j, size

        self.N = reads.shape[0]
        self.L = reads.shape[1]
        self.R = reads.shape[2]
        self.J = math.frexp(self.L)[1]-1
        for j from 0 <= j < self.J:
            size = self.L/(2**(j+1))
            self.total[j] = np.array([reads[:,k*size:(k+2)*size,:].sum(1) for k in xrange(0,2**(j+1),2)]).T
            self.valueA[j] = np.array([reads[:,k*size:(k+1)*size,:].sum(1) for k in xrange(0,2**(j+1),2)]).T
            self.valueB[j] = self.total[j] - self.valueA[j]

    def inverse_transform(self):
        """Transform a multiscale representation of the data or parameters,
        into vector representation.

        """

        if self.data:
            profile = np.array([val for k in xrange(2**self.J) \
                for val in [self.value[self.J-1][k][0],self.value[self.J-1][k][1]-self.value[self.J-1][k][0]]])
        else:
            profile = np.array([1])
            for j in xrange(self.J):
                profile = np.array([p for val in profile for p in [val,val]])
                vals = np.array([i for v in self.value[j] for i in [v,1-v]])
                profile = vals*profile

        return profile

    def copy(self):
        """ Create a copy of the class instance
        """

        cdef long j

        newcopy = Data()
        newcopy.J = self.J
        newcopy.N = self.N
        newcopy.L = self.L
        newcopy.R = self.R
        for j from 0 <= j < self.J:
            newcopy.valueA[j] = self.valueA[j]
            newcopy.valueB[j] = self.valueB[j]
            newcopy.total[j] = self.total[j]

        return newcopy


cdef class Zeta:
    """
    Inference class to store and update (E-step) the posterior
    probability that a transcription factor is bound to a motif
    instance.

    Arguments
        data : Data
        totalreads : array

    """

    def __cinit__(self, np.ndarray[np.float64_t, ndim=2] totalreads, long N, bool infer):

        cdef np.ndarray order, indices

        self.N = N
        self.total = totalreads

        if infer:
            self.prior_log_odds = np.zeros((self.N,1), dtype=float)
            self.footprint_log_likelihood_ratio = np.zeros((self.N,1), dtype=float)
            self.total_log_likelihood_ratio = np.zeros((self.N,1), dtype=float)
            self.posterior_log_odds = np.zeros((self.N,1), dtype=float)
        else:
            self.estim = np.zeros((self.N, 2),dtype=float)
            order = np.argsort(self.total.sum(1))
            indices = order[:self.N/2]
            self.estim[indices,1:] = -MAX
            indices = order[self.N/2:]
            self.estim[indices,1:] = MAX
            self.estim = np.exp(self.estim - np.max(self.estim,1).reshape(self.N,1))
            self.estim = self.estim / insum(self.estim,[1])

    cdef update(self, Data data, np.ndarray[np.float64_t, ndim=2] scores, \
        Pi pi, Tau tau, Alpha alpha, Beta beta, Omega omega, \
        Pi pi_null, Tau tau_null, str model):

        cdef long j
        cdef np.ndarray[np.float64_t, ndim=2] footprint_logodds, prior_logodds, negbin_logodds
        cdef Data lhoodA, lhoodB

        footprint_logodds = np.zeros((self.N,1), dtype=float)
        lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

        for j from 0 <= j < data.J:
            footprint_logodds += insum(lhoodA.valueA[j] - lhoodB.valueA[j],[1])

        prior_logodds = insum(beta.estim * scores, [1])
        negbin_logodds = insum(gammaln(self.total + alpha.estim.T[1]) \
                - gammaln(self.total + alpha.estim.T[0]) \
                + gammaln(alpha.estim.T[0]) - gammaln(alpha.estim.T[1]) \
                + alpha.estim.T[1] * nplog(omega.estim.T[1]) - alpha.estim.T[0] * nplog(omega.estim.T[0]) \
                + self.total * (nplog(1 - omega.estim.T[1]) - nplog(1 - omega.estim.T[0])),[1])

        self.estim[:,1:] = prior_logodds + footprint_logodds + negbin_logodds
        self.estim[:,0] = 0.
        self.estim[self.estim==np.inf] = MAX
        self.estim = np.exp(self.estim-np.max(self.estim,1).reshape(self.N,1))
        self.estim = self.estim/insum(self.estim,[1])

    cdef infer(self, Data data, np.ndarray[np.float64_t, ndim=2] scores, \
        Pi pi, Tau tau, Alpha alpha, Beta beta, Omega omega, \
        Pi pi_null, Tau tau_null, str model):

        cdef long j
        cdef Data lhoodA, lhoodB

        lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

        for j from 0 <= j < data.J:
            self.footprint_log_likelihood_ratio += insum(lhoodA.valueA[j] - lhoodB.valueA[j],[1])
        self.footprint_log_likelihood_ratio = self.footprint_log_likelihood_ratio / np.log(10)

        self.prior_log_odds = insum(beta.estim * scores, [1]) / np.log(10)

        self.total_log_likelihood_ratio = insum(gammaln(self.total + alpha.estim.T[1]) \
            - gammaln(self.total + alpha.estim.T[0]) \
            + gammaln(alpha.estim.T[0]) - gammaln(alpha.estim.T[1]) \
            + alpha.estim.T[1] * nplog(omega.estim.T[1]) - alpha.estim.T[0] * nplog(omega.estim.T[0]) \
            + self.total * (nplog(1 - omega.estim.T[1]) - nplog(1 - omega.estim.T[0])),[1])
        self.total_log_likelihood_ratio = self.total_log_likelihood_ratio / np.log(10)

        self.posterior_log_odds = self.prior_log_odds \
            + self.footprint_log_likelihood_ratio \
            + self.total_log_likelihood_ratio


cdef class Pi:
    """
    Class to store and update (M-step) the parameter `p` in the
    msCentipede model. It is also used for the parameter `p_o` in
    the msCentipede-flexbg model.

    Arguments
        J : int
        number of scales

    """

    def __cinit__(self, long J):

        cdef long j
        self.J = J
        self.value = dict()
        for j from 0 <= j < self.J:
            self.value[j] = np.empty((2**j,), dtype='float')

    def __reduce__(self):
        return (rebuild_Pi, (self.J,self.value))

    def update(self, Data data, Zeta zeta, Tau tau):
        """Update the estimates of parameter `p` (and `p_o`) in the model.
        """

        zetaestim = zeta.estim[:,1].sum()

        # call optimizer
        for j in xrange(self.J):

            # initialize optimization variable
            xo = self.value[j].copy()
            X = xo.size

            # set constraints for optimization variable
            xmin = 1./tau.estim[j]*np.ones((X,1),dtype=float)
            xmax = (1-1./tau.estim[j])*np.ones((X,1),dtype=float)
            G = np.vstack((np.diag(-1*np.ones((X,), dtype=float)), \
                    np.diag(np.ones((X,), dtype=float))))
            h = np.vstack((-1*xmin,xmax))

            # additional arguments
            args = dict([('G',G),('h',h),('data',data),('zeta',zeta),('tau',tau),('zetaestim',zetaestim),('j',j)])
    
            # call optimizer
            x_final = optimizer(xo, pi_function_gradient, pi_function_gradient_hessian, args)

            if np.isnan(x_final).any():
                print "Nan in Pi"
                raise ValueError

            if np.isinf(x_final).any():
                print "Inf in Pi"
                raise ValueError

            # store optimum in data structure
            self.value[j] = x_final

def rebuild_Pi(J, value):

    pi = Pi(J)
    pi.value = value
    return pi

cpdef tuple pi_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args):

    """Computes part of the likelihood function that has
    terms containing `pi`, along with its gradient
    """

    cdef Data data
    cdef Zeta zeta
    cdef Tau tau
    cdef long j, J, r
    cdef double f, zetaestim
    cdef np.ndarray func, F, Df, df, alpha, beta, data_alpha, data_beta

    data = args['data']
    zeta = args['zeta']
    tau = args['tau']
    zetaestim = args['zetaestim']
    j = args['j']

    J = 2**j
    func = np.zeros((data.N,J), dtype=float)
    df = np.zeros((data.N,J), dtype=float)
    alpha = x * tau.estim[j]
    beta = (1-x) * tau.estim[j]
    
    for r from 0 <= r < data.R:
        data_alpha = data.valueA[j][r] + alpha
        data_beta = data.valueB[j][r] + beta
        func += gammaln(data_alpha) + gammaln(data_beta)
        df += digamma(data_alpha) - digamma(data_beta)
    
    F = np.sum(func,1) - np.sum(gammaln(alpha) + gammaln(beta)) * data.R
    Df = tau.estim[j] * (np.sum(zeta.estim[:,1:] * df,0) \
        - zetaestim * (digamma(alpha) - digamma(beta)) * data.R)
    
    f = -1. * np.sum(zeta.estim[:,1] * F)
    Df = -1. * Df
    
    return f, Df

cpdef tuple pi_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args):

    """Computes part of the likelihood function that has
    terms containing `pi`, along with its gradient and hessian
    """

    cdef Data data
    cdef Zeta zeta
    cdef Tau tau
    cdef long j, J, r
    cdef double f, zetaestim
    cdef np.ndarray func, F, Df, df, hf, hess, Hf 

    data = args['data']
    zeta = args['zeta']
    tau = args['tau']
    zetaestim = args['zetaestim']
    j = args['j']

    hess = np.zeros((x.size,), dtype=float)

    J = 2**j
    func = np.zeros((data.N,J), dtype=float)
    df = np.zeros((data.N,J), dtype=float)
    hf = np.zeros((data.N,J), dtype=float)
    alpha = x * tau.estim[j]
    beta = (1-x) * tau.estim[j]

    for r from 0 <= r < data.R:
        data_alpha = data.valueA[j][r] + alpha
        data_beta = data.valueB[j][r] + beta
        func += gammaln(data_alpha) + gammaln(data_beta)
        df += digamma(data_alpha) - digamma(data_beta)
        hf += polygamma(1, data_alpha) + polygamma(1, data_beta)

    F = np.sum(func,1) - np.sum(gammaln(alpha) + gammaln(beta)) * data.R
    Df = tau.estim[j] * (np.sum(zeta.estim[:,1:] * df,0) \
        - zetaestim * (digamma(alpha) - digamma(beta)) * data.R)
    hess = tau.estim[j]**2 * (np.sum(zeta.estim[:,1:] * hf,0) \
        - zetaestim * (polygamma(1, alpha) + polygamma(1, beta)) * data.R)

    f = -1. * np.sum(zeta.estim[:,1] * F)
    Df = -1. * Df
    Hf = np.diag(-1.*hess)
 
    return f, Df, Hf


cdef class Tau:
    """
    Class to store and update (M-step) the parameter `tau` in the
    msCentipede model. It is also used for the parameter `tau_o` in
    the msCentipede-flexbg model.

    Arguments
        J : int
        number of scales

    """

    def __cinit__(self, long J):

        self.J = J
        self.estim = np.empty((self.J,), dtype='float')

    def __reduce__(self):
        return (rebuild_Tau, (self.J,self.estim))

    def update(self, Data data, Zeta zeta, Pi pi):
        """Update the estimates of parameter `tau` (and `tau_o`) in the model.
        """

        zetaestim = np.sum(zeta.estim[:,1])

        for j in xrange(self.J):

            # initialize optimization variables
            xo = self.estim[j:j+1]

            # set constraints for optimization variables
            minj = 1./min([np.min(pi.value[j]), np.min(1-pi.value[j])])
            xmin = np.array([minj])
            G = np.diag(-1 * np.ones((1,), dtype=float))
            h = -1*xmin.reshape(1,1)

            # additional arguments
            args = dict([('j',j),('G',G),('h',h),('data',data),('zeta',zeta),('pi',pi),('zetaestim',zetaestim)])

            # call optimizer
            try:
                x_final = optimizer(xo, tau_function_gradient, tau_function_gradient_hessian, args)
            except ValueError:
                xo = xmin+100*np.random.rand()
                bounds = [(minj, None)]
                solution = spopt.fmin_l_bfgs_b(tau_function_gradient, xo, \
                    args=(args,), bounds=bounds)
                x_final = solution[0]

            self.estim[j:j+1] = x_final

        if np.isnan(self.estim).any():
            print "Nan in Tau"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Tau"
            raise ValueError

def rebuild_Tau(J, estim):

    tau = Tau(J)
    tau.estim = estim
    return tau

cpdef tuple tau_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `tau`, and its gradient.
    """

    cdef Data data
    cdef Zeta zeta
    cdef Pi pi
    cdef long j, J, r, left, right
    cdef double F, ffunc, dff, zetaestim
    cdef np.ndarray func, ff, Df, df, alpha, beta, data_alpha, data_beta, data_x

    data = args['data']
    zeta = args['zeta']
    pi = args['pi']
    zetaestim = args['zetaestim']
    j = args['j']

    func = np.zeros((zeta.N,), dtype=float)
    ffunc = 0
    Df = np.zeros((x.size,), dtype=float)

    alpha = pi.value[j] * x
    beta = (1 - pi.value[j]) * x
    ffunc = ffunc + data.R * np.sum(gammaln(x) - gammaln(alpha) - gammaln(beta))
    dff = data.R * np.sum(digamma(x) - pi.value[j] * digamma(alpha) - (1 - pi.value[j]) * digamma(beta))
    df = np.zeros((zeta.N,), dtype=float)
    # loop over replicates
    for r from 0 <= r < data.R:

        data_alpha = data.valueA[j][r] + alpha
        data_beta = data.valueB[j][r] + beta
        data_x = data.total[j][r] + x

        func = func + np.sum(gammaln(data_alpha),1) \
            + np.sum(gammaln(data_beta),1) - np.sum(gammaln(data_x),1)

        df = df + np.sum(pi.value[j]*digamma(data_alpha),1) \
            + np.sum((1-pi.value[j])*digamma(data_beta),1) \
            - np.sum(digamma(data_x),1)

    Df[0] = -1. * (np.sum(zeta.estim[:,1] * df) + zetaestim * dff)
    F = -1. * (np.sum(zeta.estim[:,1] * func) + zetaestim * ffunc)
    
    return F, Df

cpdef tuple tau_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `tau`, and its gradient and hessian.
    """

    cdef Data data
    cdef Zeta zeta
    cdef Pi pi
    cdef long j, r, left, right
    cdef double F, ffunc, dff, hff, zetaestim
    cdef np.ndarray func, ff, Df, df, hf, hess, Hf, alpha, beta, data_alpha, data_beta, data_x

    data = args['data']
    zeta = args['zeta']
    pi = args['pi']
    zetaestim = args['zetaestim']
    j = args['j']

    func = np.zeros((zeta.N,), dtype=float)
    ffunc = 0
    Df = np.zeros((x.size,), dtype=float)
    hess = np.zeros((x.size,), dtype=float)
    # loop over each scale

    alpha = pi.value[j] * x
    beta = (1 - pi.value[j]) * x
    ffunc = ffunc + data.R * np.sum(gammaln(x) - gammaln(alpha) - gammaln(beta))
    dff = data.R * np.sum(digamma(x) - pi.value[j] * digamma(alpha) - (1 - pi.value[j]) * digamma(beta))
    hff = data.R * np.sum(polygamma(1, x) - pi.value[j]**2 * polygamma(1, alpha) \
        - (1-pi.value[j])**2 * polygamma(1, beta))
    df = np.zeros((zeta.N,), dtype=float)
    hf = np.zeros((zeta.N,), dtype=float)
    # loop over replicates
    for r from 0 <= r < data.R:

        data_alpha = data.valueA[j][r] + alpha
        data_beta = data.valueB[j][r] + beta
        data_x = data.total[j][r] + x

        func = func + np.sum(gammaln(data_alpha),1) + np.sum(gammaln(data_beta),1) - np.sum(gammaln(data_x),1)

        df = df + np.sum(pi.value[j]*digamma(data_alpha),1) \
            + np.sum((1-pi.value[j])*digamma(data_beta),1) \
            - np.sum(digamma(data_x),1)

        hf = hf + np.sum(pi.value[j]**2 * polygamma(1,data_alpha),1) \
            + np.sum((1 - pi.value[j])**2 * polygamma(1,data_beta),1) \
            - np.sum(polygamma(1,data_x),1)

    Df[0] = -1 * (np.sum(zeta.estim[:,1] * df) + zetaestim * dff)
    hess[0] = -1 * (np.sum(zeta.estim[:,1] * hf) + zetaestim * hff)
    F = -1. * (np.sum(zeta.estim[:,1] * func) + zetaestim * ffunc)
    Hf = np.diag(hess)

    return F, Df, Hf


cdef class Alpha:
    """
    Class to store and update (M-step) the parameter `alpha` in negative 
    binomial part of the msCentipede model. There is a separate parameter
    for bound and unbound states, for each replicate.

    Arguments
        R : int
        number of replicate measurements

    """

    def __cinit__(self, long R):

        self.R = R
        self.estim = np.random.rand(self.R,2)*10

    def __reduce__(self):
        return (rebuild_Alpha, (self.R,self.estim))

    def update(self, Zeta zeta, Omega omega):
        """Update the estimates of parameter `alpha` in the model.
        """

        cdef np.ndarray zetaestim, constant, xo, G, h, x_final

        zetaestim = np.sum(zeta.estim,0)
        constant = zetaestim*nplog(omega.estim)

        # initialize optimization variables
        xo = self.estim.ravel()

        # set constraints for optimization variables
        G = np.diag(-1 * np.ones((2*self.R,), dtype=float))
        h = np.zeros((2*self.R,1), dtype=float)

        args = dict([('G',G),('h',h),('omega',omega),('zeta',zeta),('constant',constant),('zetaestim',zetaestim)])

        # call optimizer
        x_final = optimizer(xo, alpha_function_gradient, alpha_function_gradient_hessian, args)
        self.estim = x_final.reshape(self.R,2)

        if np.isnan(self.estim).any():
            print "Nan in Alpha"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Alpha"
            raise ValueError

def rebuild_Alpha(R, estim):

    alpha = Alpha(R)
    alpha.estim = estim
    return alpha

cpdef tuple alpha_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `alpha`, and its gradient
    """

    cdef long r
    cdef double f, func
    cdef Zeta zeta
    cdef Omega omega
    cdef np.ndarray df, Df, constant, zetaestim, xzeta

    zeta = args['zeta']
    omega = args['omega']
    constant = args['constant']
    zetaestim = args['zetaestim']

    func = 0
    df = np.zeros((2*omega.R,), dtype='float')

    for r from 0 <= r < omega.R:
        xzeta = zeta.total[:,r:r+1] + x[2*r:2*r+2]
        func = func + np.sum(np.sum(gammaln(xzeta) * zeta.estim, 0) \
            - gammaln(x[2*r:2*r+2]) * zetaestim + constant[r] * x[2*r:2*r+2])
        df[2*r:2*r+2] = np.sum(digamma(xzeta) * zeta.estim, 0) \
            - digamma(x[2*r:2*r+2]) * zetaestim + constant[r]

    f = -1.*func
    Df = -1. * df

    return f, Df

cpdef tuple alpha_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `alpha`, and its gradient and hessian
    """

    cdef long r
    cdef double f, func
    cdef Zeta zeta
    cdef Omega omega
    cdef np.ndarray df, Df, hf, Hf, constant, zetaestim, xzeta

    zeta = args['zeta']
    omega = args['omega']
    zetaestim = args['zetaestim']
    constant = args['constant']
    
    func = 0
    df = np.zeros((2*omega.R,), dtype='float')
    hess = np.zeros((2*omega.R,), dtype='float')

    for r from 0 <= r < omega.R:
        xzeta = zeta.total[:,r:r+1] + x[2*r:2*r+2]
        func = func + np.sum(np.sum(gammaln(xzeta) * zeta.estim, 0) \
            - gammaln(x[2*r:2*r+2]) * zetaestim + constant[r] * x[2*r:2*r+2])
        df[2*r:2*r+2] = np.sum(digamma(xzeta) * zeta.estim, 0) \
            - digamma(x[2*r:2*r+2]) * zetaestim + constant[r]
        hess[2*r:2*r+2] = np.sum(polygamma(1, xzeta) * zeta.estim, 0) \
            - polygamma(1, x[2*r:2*r+2]) * zetaestim

    f = -1. * func
    Df = -1. * df
    Hf = -1. * np.diag(hess)

    return f, Df, Hf


cdef class Omega:
    """
    Class to store and update (M-step) the parameter `omega` in negative 
    binomial part of the msCentipede model. There is a separate parameter
    for bound and unbound states, for each replicate.

    Arguments
        R : int
        number of replicate measurements

    """

    def __cinit__(self, long R):

        self.R = R
        self.estim = np.random.rand(self.R,2)
        self.estim[:,1] = self.estim[:,1]/100

    def __reduce__(self):
        return (rebuild_Omega, (self.R,self.estim))

    cdef update(self, Zeta zeta, Alpha alpha):
        """Update the estimates of parameter `omega` in the model.
        """

        cdef np.ndarray numerator, denominator

        numerator = np.sum(zeta.estim,0) * alpha.estim
        denominator = np.array([np.sum(zeta.estim * (estim + zeta.total[:,r:r+1]), 0) \
            for r,estim in enumerate(alpha.estim)])
        self.estim = numerator / denominator

        if np.isnan(self.estim).any():
            print "Nan in Omega"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Omega"
            raise ValueError

def rebuild_Omega(R, estim):

    omega = Omega(R)
    omega.estim = estim
    return omega

cdef class Beta:
    """
    Class to store and update (M-step) the parameter `beta` in the logistic
    function in the prior of the msCentipede model.

    Arguments
        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

    """

    def __cinit__(self, long S):
    
        self.S = S
        self.estim = np.random.rand(self.S)

    def __reduce__(self):
        return (rebuild_Beta, (self.S,self.estim))

    def update(self, np.ndarray[np.float64_t, ndim=2] scores, Zeta zeta):
        """Update the estimates of parameter `beta` in the model.
        """

        xo = self.estim.copy()
        args = dict([('scores',scores),('zeta',zeta)])

        try:
            self.estim = optimizer(xo, beta_function_gradient, beta_function_gradient_hessian, args)
        except (ValueError, OverflowError):
            pass

        if np.isnan(self.estim).any():
            print "Nan in Beta"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Beta"
            raise ValueError

def rebuild_Beta(S, estim):

    beta = Beta(S)
    beta.estim = estim
    return beta

cpdef tuple beta_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `beta`, and its gradient.
    """

    scores = args['scores']
    zeta = args['zeta']

    arg = insum(x * scores,[1])
    
    func = arg * zeta.estim[:,1:] - nplog(1 + np.exp(arg))
    f = -1. * func.sum()
    
    Df = -1 * np.sum(scores * (zeta.estim[:,1:] - logistic(-arg)),0)
    
    return f, Df

cpdef tuple beta_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args):
    """Computes part of the likelihood function that has
    terms containing `beta`, and its gradient and hessian.
    """

    scores = args['scores']
    zeta = args['zeta']

    arg = insum(x * scores,[1])
    
    func = arg * zeta.estim[:,1:] - nplog(1 + np.exp(arg))
    f = -1. * func.sum()
    
    Df = -1 * np.sum(scores * (zeta.estim[:,1:] - logistic(-arg)),0)
    
    larg = scores * logistic(arg) * logistic(-arg)
    Hf = np.dot(scores.T, larg)
    
    return f, Df, Hf


def optimizer(np.ndarray[np.float64_t, ndim=1] xo, function_gradient, function_gradient_hessian, dict args):
    """Calls the appropriate nonlinear convex optimization solver 
    in the package `cvxopt` to find optimal values for the relevant
    parameters, given subroutines that evaluate a function, 
    its gradient, and hessian, this subroutine 

    Arguments
        function : function object
        evaluates the function at the specified parameter values

        gradient : function object
        evaluates the gradient of the function

        hessian : function object
        evaluates the hessian of the function

    """

    def F(x=None, z=None):
        """A subroutine that the cvxopt package can call to get 
        values of the function, gradient and hessian during
        optimization.
        """

        if x is None:
            return 0, cvx.matrix(x_init)

        xx = np.array(x).ravel().astype(np.float64)

        if z is None:

            # compute likelihood function and gradient
            f, Df = function_gradient(xx, args)

            # check for infs and nans in function and gradient
            if np.isnan(f) or np.isinf(f):
                f = np.array([np.finfo('float32').max]).astype('float')
            else:
                f = np.array([f]).astype('float')
            if np.isnan(Df).any() or np.isinf(Df).any():
                Df = -1 * np.finfo('float32').max * np.ones((1,xx.size), dtype=float)
            else:
                Df = Df.reshape(1,xx.size)

            return cvx.matrix(f), cvx.matrix(Df)

        else:

            # compute likelihood function, gradient, and hessian
            f, Df, hess = function_gradient_hessian(xx, args)

            # check for infs and nans in function and gradient
            if np.isnan(f) or np.isinf(f):
                f = np.array([np.finfo('float32').max]).astype('float')
            else:
                f = np.array([f]).astype('float')
            if np.isnan(Df).any() or np.isinf(Df).any():
                Df = -1 * np.finfo('float32').max * np.ones((1,xx.size), dtype=float)
            else:
                Df = Df.reshape(1,xx.size)

            Hf = z[0] * hess
            return cvx.matrix(f), cvx.matrix(Df), cvx.matrix(Hf)

    # warm start for the optimization
    V = xo.size
    x_init = xo.reshape(V,1)

    # call the optimization subroutine in cvxopt
    if args.has_key('G'):
        # call a constrained nonlinear solver
        solution = solvers.cp(F, G=cvx.matrix(args['G']), h=cvx.matrix(args['h']))
    else:
        # call an unconstrained nonlinear solver
        solution = solvers.cp(F)

    x_final = np.array(solution['x']).ravel()

    return x_final


cdef tuple compute_footprint_likelihood(Data data, Pi pi, Tau tau, Pi pi_null, Tau tau_null, str model):
    """Evaluates the likelihood function for the 
    footprint part of the bound model and background model.

    Arguments
        data : Data
        transformed read count data 

        pi : Pi
        estimate of mean footprint parameters at bound sites

        tau : Tau
        estimate of footprint heterogeneity at bound sites

        pi_null : Pi
        estimate of mean cleavage pattern at unbound sites

        tau_null : Tau or None
        estimate of cleavage heterogeneity at unbound sites

        model : string 
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

    """

    cdef long j, r
    cdef np.ndarray valueA, valueB
    cdef Data lhood_bound, lhood_unbound

    lhood_bound = Data()
    lhood_unbound = Data()

    for j from 0 <= j < data.J:
        valueA = np.sum(data.valueA[j],0)
        valueB = np.sum(data.valueB[j],0)
        
        lhood_bound.valueA[j] = np.sum([gammaln(data.valueA[j][r] + pi.value[j] * tau.estim[j]) \
            + gammaln(data.valueB[j][r] + (1 - pi.value[j]) * tau.estim[j]) \
            - gammaln(data.total[j][r] + tau.estim[j]) + gammaln(tau.estim[j]) \
            - gammaln(pi.value[j] * tau.estim[j]) - gammaln((1 - pi.value[j]) * tau.estim[j]) \
            for r in xrange(data.R)],0)

        if model in ['msCentipede','msCentipede_flexbgmean']:
            
            lhood_unbound.valueA[j] = valueA * nplog(pi_null.value[j]) \
                + valueB * nplog(1 - pi_null.value[j])
    
        elif model=='msCentipede_flexbg':
            
            lhood_unbound.valueA[j] = np.sum([gammaln(data.valueA[j][r] + pi_null.value[j] * tau_null.estim[j]) \
                + gammaln(data.valueB[j][r] + (1 - pi_null.value[j]) * tau_null.estim[j]) \
                - gammaln(data.total[j][r] + tau_null.estim[j]) + gammaln(tau_null.estim[j]) \
                - gammaln(pi_null.value[j] * tau_null.estim[j]) - gammaln((1 - pi_null.value[j]) * tau_null.estim[j]) \
                for r in xrange(data.R)],0)

    return lhood_bound, lhood_unbound


cdef double likelihood(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
    Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
    Omega omega, Pi pi_null, Tau tau_null, str model):
    """Evaluates the likelihood function of the full
    model, given estimates of model parameters.

    Arguments
        data : Data
        transformed read count data

        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

        zeta : zeta
        expected value of factor binding state for each site.

        pi : Pi
        estimate of mean footprint parameters at bound sites

        tau : Tau
        estimate of footprint heterogeneity at bound sites

        alpha : Alpha
        estimate of negative binomial parameters for each replicate

        beta : Beta
        weights for various scores in the logistic function 

        omega : Omega
        estimate of negative binomial parameters for each replicate

        pi_null : Pi
        estimate of mean cleavage pattern at unbound sites

        tau_null : Tau or None
        estimate of cleavage heterogeneity at unbound sites

        model : string
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

    """

    cdef long j
    cdef double L
    cdef np.ndarray apriori, footprint, null, P_1, P_0, LL
    cdef Data lhoodA, lhoodB

    apriori = insum(beta.estim * scores,[1])

    lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

    footprint = np.zeros((data.N,1),dtype=float)
    for j from 0 <= j < data.J:
        footprint += insum(lhoodA.valueA[j],[1])

    P_1 = footprint + insum(gammaln(zeta.total + alpha.estim[:,1]) - gammaln(alpha.estim[:,1]) \
        + alpha.estim[:,1] * nplog(omega.estim[:,1]) + zeta.total * nplog(1 - omega.estim[:,1]), [1])
    P_1[P_1==np.inf] = MAX
    P_1[P_1==-np.inf] = -MAX

    null = np.zeros((data.N,1), dtype=float)
    for j from 0 <= j < data.J:
        null += insum(lhoodB.valueA[j],[1])

    P_0 = null + insum(gammaln(zeta.total + alpha.estim[:,0]) - gammaln(alpha.estim[:,0]) \
        + alpha.estim[:,0] * nplog(omega.estim[:,0]) + zeta.total * nplog(1 - omega.estim[:,0]), [1])
    P_0[P_0==np.inf] = MAX
    P_0[P_0==-np.inf] = -MAX

    LL = P_0 * zeta.estim[:,:1] + P_1 * zeta.estim[:,1:] + apriori * (1 - zeta.estim[:,:1]) \
        - nplog(1 + np.exp(apriori)) - insum(zeta.estim * nplog(zeta.estim),[1])
 
    L = LL.sum() / data.N

    if np.isnan(L):
        print "Nan in LogLike"
        return -np.inf

    if np.isinf(L):
        print "Inf in LogLike"
        return -np.inf

    return L


cdef EM(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
    Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
    Omega omega, Pi pi_null, Tau tau_null, str model):
    """This subroutine updates all model parameters once and computes an
    estimate of the posterior probability of binding.

    Arguments
        data : Data
        transformed read count data

        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

        zeta : zeta
        expected value of factor binding state for each site.

        pi : Pi
        estimate of mean footprint parameters at bound sites

        tau : Tau
        estimate of footprint heterogeneity at bound sites

        alpha : Alpha
        estimate of negative binomial parameters for each replicate

        beta : Beta
        weights for various scores in the logistic function 

        omega : Omega
        estimate of negative binomial parameters for each replicate

        pi_null : Pi
        estimate of mean cleavage pattern at unbound sites

        tau_null : Tau or None
        estimate of cleavage heterogeneity at unbound sites

        model : string
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

    """

    cdef double starttime

    # update binding posteriors
    zeta.update(data, scores, pi, tau, \
            alpha, beta, omega, pi_null, tau_null, model)

    # update multi-scale parameters
    #starttime = time.time()
    pi.update(data, zeta, tau)
    #print "p_jk update in %.3f secs"%(time.time()-starttime)

    #starttime = time.time()
    tau.update(data, zeta, pi)
    #print "tau update in %.3f secs"%(time.time()-starttime)

    # update negative binomial parameters
    #starttime = time.time()
    omega.update(zeta, alpha)
    #print "omega update in %.3f secs"%(time.time()-starttime)

    #starttime = time.time()
    alpha.update(zeta, omega)
    #print "alpha update in %.3f secs"%(time.time()-starttime)

    # update prior parameters
    #starttime = time.time()
    beta.update(scores, zeta)
    #print "beta update in %.3f secs"%(time.time()-starttime)


cdef square_EM(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
    Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
    Omega omega, Pi pi_null, Tau tau_null, str model):
    """Accelerated update of model parameters and posterior probability of binding.

    Arguments
        data : Data
        transformed read count data

        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

        zeta : zeta
        expected value of factor binding state for each site.

        pi : Pi
        estimate of mean footprint parameters at bound sites

        tau : Tau
        estimate of footprint heterogeneity at bound sites

        alpha : Alpha
        estimate of negative binomial parameters for each replicate

        beta : Beta
        weights for various scores in the logistic function 

        omega : Omega
        estimate of negative binomial parameters for each replicate

        pi_null : Pi
        estimate of mean cleavage pattern at unbound sites

        tau_null : Tau or None
        estimate of cleavage heterogeneity at unbound sites

        model : string
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

    """

    cdef long j, step
    cdef double a
    cdef np.ndarray r, v, varA, varB, varC, invalid, newparam
    cdef list parameters, oldvar, oldvars, R, V

    parameters = [pi, tau, alpha, omega]
    oldvar = []
    for parameter in parameters:
        try:
            oldvar.append(parameter.estim.copy())
        except AttributeError:
            oldvar.append(np.hstack([parameter.value[j].copy() for j in xrange(parameter.J)]))
    oldvars = [oldvar]

    # take two update steps
    for step in [0,1]:
        EM(data, scores, zeta, pi, tau, alpha, beta, omega, pi_null, tau_null, model)
        oldvar = []
        for parameter in parameters:
            try:
                oldvar.append(parameter.estim.copy())
            except AttributeError:
                oldvar.append(np.hstack([parameter.value[j].copy() for j in xrange(parameter.J)]))
        oldvars.append(oldvar)

    R = [oldvars[1][j]-oldvars[0][j] for j in xrange(len(parameters))]
    V = [oldvars[2][j]-oldvars[1][j]-R[j] for j in xrange(len(parameters))]
    a = -1.*np.sqrt(np.sum([(r*r).sum() for r in R]) / np.sum([(v*v).sum() for v in V]))

    if a>-1:
        a = -1.

    # given two update steps, compute an optimal step that achieves
    # a better likelihood than the two steps.
    a_ok = False
    while not a_ok:
        invalid = np.zeros((0,), dtype='bool')
        for parameter,varA,varB,varC in zip(parameters,oldvars[0],oldvars[1],oldvars[2]):
            try:
                parameter.estim = (1+a)**2*varA - 2*a*(1+a)*varB + a**2*varC
                # ensure constraints on variables are satisfied
                invalid = np.hstack((invalid,(parameter.estim<=0).ravel()))
            except AttributeError:
                newparam = (1+a)**2*varA - 2*a*(1+a)*varB + a**2*varC
                # ensure constraints on variables are satisfied
                invalid = np.hstack((invalid, np.logical_or(newparam<0, newparam>1)))
                parameter.value = dict([(j,newparam[2**j-1:2**(j+1)-1]) \
                    for j in xrange(parameter.J)])
        if np.any(invalid):
            a = (a-1)/2.
            if np.abs(a+1)<1e-4:
                a = -1.
        else:
            a_ok = True

    EM(data, scores, zeta, pi, tau, alpha, beta, omega, pi_null, tau_null, model)


def estimate_optimal_model(np.ndarray[np.float64_t, ndim=3] reads, \
    np.ndarray[np.float64_t, ndim=2] totalreads, \
    np.ndarray[np.float64_t, ndim=2] scores, \
    np.ndarray[np.float64_t, ndim=3] background, \
    str model, str log_file, long restarts, double mintol):
    """Learn the model parameters by running an EM algorithm till convergence.
    Return the optimal parameter estimates from a number of EM results starting 
    from random restarts.

    Arguments
        reads : array
        array of read counts at each base in a genomic window,
        across motif instances and several measurement replicates.

        totalreads : array
        array of total read counts in a genomic window,
        across motif instances and several measurement replicates.
        the size of the genomic window can be different for 
        `reads` and `totalreads`.

        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

        background : array
        a uniform, normalized array for a uniform background model.
        when sequencing reads from genomic DNA are available, this
        is an array of read counts at each base in a genomic window,
        across motif instances.

        model : string
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

        restarts : int
        number of independent runs of model learning

        mintol : float
        convergence criterion

    """

    cdef long restart, iteration, err
    cdef double change, maxLoglike, Loglike, tol, itertime, totaltime
    cdef np.ndarray oldtau, negbinmeans
    cdef Data data, data_null
    cdef Beta beta
    cdef Alpha alpha
    cdef Omega omega
    cdef Pi pi, pi_null
    cdef Tau tau, tau_null
    cdef Zeta zeta, zeta_null

    log = "transforming data into multiscale representation ..."
    log_handle = open(log_file,'a')
    log_handle.write(log)
    log_handle.close()
    print log

    # transform data into multiscale representation
    data = Data()
    data.transform_to_multiscale(reads)
    data_null = Data()
    data_null.transform_to_multiscale(background)
    del reads

    # transform matrix of PWM scores and other prior information
    scores = np.hstack((np.ones((data.N,1), dtype=float), scores))
    S = scores.shape[1]

    log_handle = open(log_file,'a')
    log_handle.write("done\n")
    log_handle.close()

    # set background model
    pi_null = Pi(data_null.J)
    for j in xrange(pi_null.J):
        pi_null.value[j] = np.sum(np.sum(data_null.valueA[j],0),0) / np.sum(np.sum(data_null.total[j],0),0).astype('float')
    
    tau_null = Tau(data_null.J)
    if model=='msCentipede_flexbg':

        log = "learning a flexible background model ..."
        log_handle = open(log_file,'a')
        log_handle.write(log)
        log_handle.close()
        print log

        zeta_null = Zeta(background.sum(1), data_null.N, False)
        zeta_null.estim[:,1] = 1
        zeta_null.estim[:,0] = 0

        # iterative update of background model; 
        # evaluate convergence based on change in estimated
        # background overdispersion
        change = np.inf
        while change>1e-2:
            oldtau = tau_null.estim.copy()
            
            tau_null.update(data_null, zeta_null, pi_null)
            pi_null.update(data_null, zeta_null, tau_null)

            change = np.abs(oldtau-tau_null.estim).sum() / tau_null.J

        log_handle = open(log_file,'a')
        log_handle.write("done\n")
        log_handle.close()

    maxLoglike = -np.inf
    restart = 0
    err = 1
    while restart<restarts:

        totaltime = time.time()
        log = "starting model estimation (restart %d)"%(restart+1)
        log_handle = open(log_file,'a')
        log_handle.write(log+'\n')
        log_handle.close()
        print log

        # initialize multi-scale model parameters
        pi = Pi(data.J)
        tau = Tau(data.J)

        # initialize negative binomial parameters
        alpha = Alpha(data.R)
        omega = Omega(data.R)

        # initialize prior parameters
        beta = Beta(S)

        # initialize posterior over latent variables
        zeta = Zeta(totalreads, data.N, False)
        for j in xrange(pi.J):
            pi.value[j] = np.sum(data.valueA[j][0] * zeta.estim[:,1:],0) \
                / np.sum(data.total[j][0] * zeta.estim[:,1:],0).astype('float')
            mask = pi.value[j]>0
            pi.value[j][~mask] = pi.value[j][mask].min()
            mask = pi.value[j]<1
            pi.value[j][~mask] = pi.value[j][mask].max()
            minj = 1./min([pi.value[j].min(), (1-pi.value[j]).min()])
            if minj<2:
                minj = 2.
            tau.estim[j] = minj+10*np.random.rand()

        # initial log likelihood of the model
        Loglike = likelihood(data, scores, zeta, pi, tau, \
                alpha, beta, omega, pi_null, tau_null, model)

        log = "initial log likelihood = %.2e"%Loglike
        log_handle = open(log_file,'a')
        log_handle.write(log+'\n')
        log_handle.close()
        print log

        tol = np.inf
        iteration = 0

        while np.abs(tol)>mintol:

            itertime = time.time()
            square_EM(data, scores, zeta, pi, tau, \
                    alpha, beta, omega, pi_null, tau_null, model)

            newLoglike = likelihood(data, scores, zeta, pi, tau, \
                    alpha, beta, omega, pi_null, tau_null, model)

            tol = newLoglike - Loglike
            Loglike = newLoglike
            log = "iteration %d: log likelihood = %.2e, change in log likelihood = %.2e, iteration time = %.3f secs"%(iteration+1, Loglike, tol, time.time()-itertime)
            log_handle = open(log_file,'a')
            log_handle.write(log+'\n')
            log_handle.close()
            print log
            iteration += 1
        totaltime = (time.time()-totaltime)/60.

        # test if mean cleavage rate at bound sites is greater than at 
        # unbound sites, for each replicate; avoids local optima issues.
        negbinmeans = alpha.estim * (1-omega.estim)/omega.estim
        if np.any(negbinmeans[:,0]<negbinmeans[:,1]):
            restart += 1
            # choose these parameter estimates, if the likelihood is greater.
            if Loglike>maxLoglike:
                maxLoglikeres = Loglike
                if model in ['msCentipede','msCentipede_flexbgmean']:
                    footprint_model = (pi, tau, pi_null)
                elif model=='msCentipede_flexbg':
                    footprint_model = (pi, tau, pi_null, tau_null)
                count_model = (alpha, omega)
                prior = beta

    return footprint_model, count_model, prior


def infer_binding_posterior(reads, totalreads, scores, background, footprint, negbinparams, prior, model):
    """Infer posterior probability of factor binding, given optimal model parameters.

    Arguments
        reads : array
        array of read counts at each base in a genomic window,
        across motif instances and several measurement replicates.

        totalreads : array
        array of total read counts in a genomic window,
        across motif instances and several measurement replicates.
        the size of the genomic window can be different for 
        `reads` and `totalreads`.

        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

        background : array
        a uniform, normalized array for a uniform background model.
        when sequencing reads from genomic DNA are available, this
        is an array of read counts at each base in a genomic window,
        across motif instances.

        footprint : tuple
        (Pi, Tau) instances
        estimate of footprint model parameters

        negbinparams : tuple
        (Alpha, Omega) instances
        estimate of negative binomial model parameters

        prior : Beta
        estimate of weights in logistic function in the prior

        model : string
        {msCentipede, msCentipede-flexbgmean, msCentipede-flexbg}

    """

    data = Data()
    data.transform_to_multiscale(reads)
    data_null = Data()
    data_null.transform_to_multiscale(background)
    scores = np.hstack((np.ones((data.N,1), dtype=float), scores))
    del reads

    # negative binomial parameters
    alpha = negbinparams[0]
    omega = negbinparams[1]
    
    # weights in logistic function in the prior
    beta = prior

    # multiscale parameters
    pi = footprint[0]
    tau = footprint[1]
    
    # setting background model
    pi_null = footprint[2]
    for j in xrange(pi_null.J):
        pi_null.value[j] = np.sum(np.sum(data_null.valueA[j],0),0) \
            / np.sum(np.sum(data_null.total[j],0),0).astype('float')
    tau_null = None

    if model=='msCentipede_flexbg':

        tau_null = footprint[3]

        if data_null.N>1000:

            zeta_null = Zeta(background.sum(1), data_null.N, False)
            zeta_null.estim[:,1] = 1
            zeta_null.estim[:,0] = 0

            # iterative update of background model, when
            # accounting for overdispersion
            change = np.inf
            while change>1e-1:
                change = tau_null.estim.copy()
                
                pi_null.update(data_null, zeta_null, tau_null)

                tau_null.update(data_null, zeta_null, pi_null)

                change = np.abs(change-tau_null.estim).sum()

    zeta = Zeta(totalreads, data.N, True)

    zeta.infer(data, scores, pi, tau, alpha, beta, omega, \
        pi_null, tau_null, model)
    
    return zeta.posterior_log_odds, \
        zeta.prior_log_odds, zeta.footprint_log_likelihood_ratio, \
        zeta.total_log_likelihood_ratio

