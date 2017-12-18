import numpy as np
import scipy.optimize as spopt
import cvxopt as cvx
from cvxopt import solvers
from scipy.special import digamma, gammaln, polygamma
import time, math, pdb

# suppress optimizer output
solvers.options['show_progress'] = False
solvers.options['maxiters'] = 40
np.random.seed(10)

# defining some constants
EPS = np.finfo(np.double).tiny
MAX = np.finfo(np.double).max

# defining some simple functions
logistic = lambda x: 1./(1+np.exp(x))
insum = lambda x,axes: np.apply_over_axes(np.sum,x,axes)

def outsum(arr):
    """Summation over the first axis, without changing length of shape.

    Arguments
        arr : array

    Returns
        thesum : array

    .. note::
        This implementation is much faster than `numpy.sum`.

    """

    thesum = sum([a for a in arr])
    shape = [1]
    shape.extend(list(thesum.shape))
    thesum = thesum.reshape(tuple(shape))
    return thesum

def nplog(x):
    """Compute the natural logarithm, handling very
    small floats appropriately.

    """
    try:
        x[x<EPS] = EPS
    except TypeError:
        x = max([x,EPS])
    return np.log(x)


class Data:
    """
    A data structure to store a multiscale representation of
    chromatin accessibility read counts across `N` genomic windows of
    length `L` in `R` replicates.

    Arguments
        reads : array

    """

    def __init__(self, reads=None):

        if reads is None:
            self.N = 0
            self.L = 0
            self.R = 0
            self.J = 0
            self.value = dict()
            self.total = dict()
        else:
            self.N, self.L, self.R = reads.shape
            self.J = math.frexp(self.L)[1]-1
            self.value = dict()
            self.total = dict()
            self.transform(reads)

    def transform(self, profile):
        """Transform a vector of read counts or parameter values
        into a multiscale representation.
        
        .. note::
            See msCentipede manual for more details.   

        """

        for j in xrange(self.J):
            size = self.L/(2**(j+1))
            self.total[j] = np.array([profile[:,k*size:(k+2)*size,:].sum(1) for k in xrange(0,2**(j+1),2)]).T
            self.value[j] = np.array([profile[:,k*size:(k+1)*size,:].sum(1) for k in xrange(0,2**(j+1),2)]).T

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

        newcopy = Data()
        newcopy.J = self.J
        newcopy.N = self.N
        newcopy.L = self.L
        newcopy.R = self.R
        for j in xrange(self.J):
            newcopy.value[j] = self.value[j]
            newcopy.total[j] = self.total[j]

        return newcopy


class Zeta():
    """
    Inference class to store and update (E-step) the posterior
    probability that a transcription factor is bound to a motif
    instance.

    Arguments
        data : Data
        totalreads : array

    """

    def __init__(self, data, totalreads, infer=False):

        self.N = data.N
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

    def update(self, data, scores, pi, tau, alpha, beta, omega, \
        pi_null, tau_null, model):

        footprint_logodds = np.zeros((self.N,1),dtype=float)
        lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

        for j in xrange(data.J):
            footprint_logodds += insum(lhoodA.value[j] - lhoodB.value[j],[1])

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

    def infer(self, data, scores, pi, tau, alpha, beta, omega, \
        pi_null, tau_null, model):

        lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

        for j in xrange(data.J):
            self.footprint_log_likelihood_ratio += insum(lhoodA.value[j] - lhoodB.value[j],[1])
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


class Pi(Data):
    """
    Class to store and update (M-step) the parameter `p` in the
    msCentipede model. It is also used for the parameter `p_o` in
    the msCentipede-flexbg model.

    Arguments
        J : int
        number of scales

    """

    def __init__(self, J):

        Data.__init__(self)
        self.J = J
        for j in xrange(self.J):
            self.value[j] = np.empty((2**j,), dtype='float')

    def update(self, data, zeta, tau):
        """Update the estimates of parameter `p` (and `p_o`) in the model.
        """

        def function(x, kwargs):
            """Computes part of the likelihood function that has
            terms containing `pi`.
            """
            data = kwargs['data']
            zeta = kwargs['zeta']
            tau = kwargs['tau']
            j = kwargs['j']

            func = np.zeros(data.value[j][0].shape, dtype=float)
            for r in xrange(data.R):
                func += gammaln(data.value[j][r] + tau.estim[j] * x) \
                    + gammaln(data.total[j][r] - data.value[j][r] + tau.estim[j] * (1-x)) \
                    - gammaln(tau.estim[j] * x) - gammaln(tau.estim[j] * (1-x))
            f = -1. * np.sum(zeta.estim[:,1] * np.sum(func,1))
            return f

        def gradient(x, kwargs):
            """Computes gradient of the likelihood function with respect to `pi`.
            """

            data = kwargs['data']
            zeta = kwargs['zeta']
            tau = kwargs['tau']
            j = kwargs['j']

            df = np.zeros(data.value[j][0].shape, dtype=float)
            for r in xrange(data.R):
                df += digamma(data.value[j][r] + tau.estim[j] * x) \
                - digamma(data.total[j][r] - data.value[j][r] + tau.estim[j] * (1-x)) \
                - digamma(tau.estim[j] * x) + digamma(tau.estim[j] * (1-x))
            Df = -1. * tau.estim[j] * np.sum(zeta.estim[:,1:] * df,0)
            return Df

        def hessian(x, kwargs):
            """Computes hessian of the likelihood function with respect to `pi`.
            """

            data = kwargs['data']
            zeta = kwargs['zeta']
            tau = kwargs['tau']
            j = kwargs['j']

            hf = np.zeros(data.value[j][0].shape, dtype=float)
            for r in xrange(data.R):
                hf += polygamma(1, data.value[j][r] + tau.estim[j] * x) \
                + polygamma(1, data.total[j][r] - data.value[j][r] + tau.estim[j] * (1-x)) \
                - polygamma(1, tau.estim[j] * x) - polygamma(1, tau.estim[j] * (1-x))
            hess = -1. * tau.estim[j]**2 * np.sum(zeta.estim[:,1:] * hf,0)
            
            Hf = np.diag(hess)
            return Hf

        for j in xrange(self.J):

            # initialize optimization variable
            xo = self.value[j]
            X = xo.size

            # set constraints for optimization variable
            if tau.estim[j]<2:
                xmin = 0.5*np.ones((X,1), dtype='float')
                xmax = 0.5*np.ones((X,1), dtype='float')
            else:
                xmin = 1./tau.estim[j]*np.ones((X,1), dtype='float')
                xmax = (1-1./tau.estim[j])*np.ones((X,1), dtype='float')
            G = np.vstack((np.diag(-1*np.ones((X,), dtype='float')), \
                    np.diag(np.ones((X,), dtype='float'))))
            h = np.vstack((-1*xmin,xmax))

            args = dict([('G',G),('h',h),('data',data),('zeta',zeta),('tau',tau),('j',j)])

            # call optimizer
            optimized = False
            while not optimized:
                try:
                    self.value[j] = optimizer(xo, function, gradient, hessian, args)
                    optimized = True
                except ValueError:
                    dx = xmax-xmin
                    xo[dx>0] = xmin + np.random.rand(X,1)/dx
                    xo[dx==0] = xmin

            if np.isnan(self.value[j]).any():
                print "Nan in Pi"
                raise ValueError

            if np.isinf(self.value[j]).any():
                print "Inf in Pi"
                raise ValueError

    def avoid_edges(self):

        for j in xrange(self.J):
            self.value[j][self.value[j]<1e-10] = 1e-10
            self.value[j][self.value[j]>1-1e-10] = 1-1e-10

class Tau():
    """
    Class to store and update (M-step) the parameter `tau` in the
    msCentipede model. It is also used for the parameter `tau_o` in
    the msCentipede-flexbg model.

    Arguments
        J : int
        number of scales

    """

    def __init__(self, J):

        self.J = J
        self.estim = np.empty((self.J,), dtype='float')

    def update(self, data, zeta, pi):
        """Update the estimates of parameter `tau` (and `tau_o`) in the model.
        """

        def function(x, kwargs):
            """Computes part of the likelihood function that has
            terms containing `tau`.
            """

            data = kwargs['data']
            zeta = kwargs['zeta']
            pi = kwargs['pi']
            j = kwargs['j']

            func = np.zeros(zeta.estim[:,1].shape, dtype=float)
            # loop over replicates
            for r in xrange(data.R):
                F = gammaln(data.value[j][r] + pi.value[j] * x) \
                    + gammaln(data.total[j][r] - data.value[j][r] + (1 - pi.value[j]) * x) \
                    - gammaln(data.total[j][r] + x) + gammaln(x) \
                    - gammaln(pi.value[j] * x) - gammaln((1 - pi.value[j]) * x)
                func += np.sum(F, 1)

            F = -1. * np.sum(zeta.estim[:,1] * func)
            return F

        def gradient(x, kwargs):
            """Computes gradient of the likelihood function with respect to `tau`.
            """

            data = kwargs['data']
            zeta = kwargs['zeta']
            pi = kwargs['pi']
            j = kwargs['j']

            # loop over replicates
            Df = np.empty((1,), dtype='float')
            df = np.zeros(zeta.estim[:,1].shape, dtype=float)
            for r in xrange(data.R):
                f = pi.value[j] * digamma(data.value[j][r] + pi.value[j] * x) \
                    + (1 - pi.value[j]) * digamma(data.total[j][r] - data.value[j][r] + (1 - pi.value[j]) * x) \
                    - digamma(data.total[j][r] + x) + digamma(x) \
                    - pi.value[j] * digamma(pi.value[j] * x) - (1 - pi.value[j]) * digamma((1 - pi.value[j]) * x)
                df += np.sum(f, 1)
            Df[0] = -1 * np.sum(zeta.estim[:,1] * df)
            return Df

        def hessian(x, kwargs):
            """Computes hessian of the likelihood function with respect to `tau`.
            """

            data = kwargs['data']
            zeta = kwargs['zeta']
            pi = kwargs['pi']
            j = kwargs['j']

            # loop over replicates
            hess = np.empty((1,), dtype='float')
            hf = np.zeros(zeta.estim[:,1].shape, dtype=float)
            for r in xrange(data.R):
                f = pi.value[j]**2 * polygamma(1, data.value[j][r] + pi.value[j] * x) \
                    + (1 - pi.value[j])**2 * polygamma(1, data.total[j][r] - data.value[j][r] + (1 - pi.value[j]) * x) \
                    - polygamma(1, data.total[j][r] + x) + polygamma(1, x) \
                    - pi.value[j]**2 * polygamma(1, pi.value[j] * x) \
                    - (1 - pi.value[j])**2 * polygamma(1, (1 - pi.value[j]) * x)
                hf += np.sum(f, 1)
            hess[0] = -1 * np.sum(zeta.estim[:,1] * hf)

            Hf = np.diag(hess)
            return Hf

        for j in xrange(self.J):

            # initialize optimization variables
            xo = self.estim[j:j+1]

            # set constraints for optimization variables
            G = np.diag(-1 * np.ones((1,), dtype=float))
            minj = 1./min([pi.value[j].min(), (1-pi.value[j]).min()])
            xmin = np.array(minj).reshape(1,1)
            h = -1*xmin

            args = dict([('j',j),('G',G),('h',h),('data',data),('zeta',zeta),('pi',pi)])

            # call optimizer
            optimized = False
            while not optimized:
                try:
                    x_final = optimizer(xo, function, gradient, hessian, args)
                    optimized = True
                except ValueError as err:
                    xo = xmin.ravel()+100*np.random.rand()
                    bounds = [(minj, None)]
                    solution = spopt.fmin_l_bfgs_b(function, xo, fprime=gradient, \
                        args=(args,), bounds=bounds)
                    x_final = solution[0]
                    optimized = True
            self.estim[j:j+1] = x_final

        if np.isnan(self.estim).any():
            print "Nan in Tau"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Tau"
            raise ValueError


class Alpha():
    """
    Class to store and update (M-step) the parameter `alpha` in negative 
    binomial part of the msCentipede model. There is a separate parameter
    for bound and unbound states, for each replicate.

    Arguments
        R : int
        number of replicate measurements

    """

    def __init__(self, R):

        self.R = R
        self.estim = np.random.rand(self.R,2)*10

    def update(self, zeta, omega):
        """Update the estimates of parameter `alpha` in the model.
        """

        def function(x, kwargs):
            """Computes part of the likelihood function that has
            terms containing `alpha`.
            """

            zeta = kwargs['zeta']
            omega = kwargs['omega']
            constant = kwargs['constant']
            zetaestim = kwargs['zetaestim']

            func = np.array([outsum(gammaln(zeta.total[:,r:r+1] + x[2*r:2*r+2]) * zeta.estim) \
                    - gammaln(x[2*r:2*r+2]) * zetaestim[0] + constant[r] * x[2*r:2*r+2] \
                    for r in xrange(omega.R)])
            f = -1.*func.sum()
            return f

        def gradient(x, kwargs):
            """Computes gradient of the likelihood function with 
            respect to `omega`.
            """

            zeta = kwargs['zeta']
            omega = kwargs['omega']
            zetaestim = kwargs['zetaestim']
            constant = kwargs['constant']

            df = []
            for r in xrange(omega.R):
                df.append(outsum(digamma(zeta.total[:,r:r+1] + x[2*r:2*r+2]) * zeta.estim)[0] \
                    - digamma(x[2*r:2*r+2]) * zetaestim[0] + constant[r])
            Df = -1. * np.hstack(df)
            return Df

        def hessian(x, kwargs):
            """Computes hessian of the likelihood function with 
            respect to `omega`.
            """

            zeta = kwargs['zeta']
            omega = kwargs['omega']
            zetaestim = kwargs['zetaestim']
            constant = kwargs['constant']
            
            hess = []
            for r in xrange(omega.R):
                hess.append(outsum(polygamma(1, zeta.total[:,r:r+1] + x[2*r:2*r+2]) * zeta.estim)[0] \
                    - polygamma(1, x[2*r:2*r+2]) * zetaestim[0])
            Hf = -1. * np.diag(np.hstack(hess))
            return Hf

        constant = [nplog(omega.estim[r]) * outsum(zeta.estim)[0] for r in xrange(self.R)]
        zetaestim = outsum(zeta.estim)

        # initialize optimization variables
        xo = self.estim.ravel()

        # set constraints for optimization variables
        G = np.diag(-1 * np.ones(xo.shape, dtype=float))
        h = np.zeros((xo.size,1), dtype=float)

        args = dict([('G',G),('h',h),('omega',omega),('zeta',zeta),('constant',constant),('zetaestim',zetaestim)])

        # call optimizer
        x_final = optimizer(xo, function, gradient, hessian, args)
        self.estim = x_final.reshape(self.estim.shape)

        if np.isnan(self.estim).any():
            print "Nan in Alpha"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Alpha"
            raise ValueError


class Omega():
    """
    Class to store and update (M-step) the parameter `omega` in negative 
    binomial part of the msCentipede model. There is a separate parameter
    for bound and unbound states, for each replicate.

    Arguments
        R : int
        number of replicate measurements

    """

    def __init__(self, R):

        self.R = R
        self.estim = np.random.rand(self.R,2)
        self.estim[:,1] = self.estim[:,1]/100

    def update(self, zeta, alpha):
        """Update the estimates of parameter `omega` in the model.
        """

        numerator = outsum(zeta.estim)[0] * alpha.estim
        denominator = np.array([outsum(zeta.estim * (estim + zeta.total[:,r:r+1]))[0] \
            for r,estim in enumerate(alpha.estim)])
        self.estim = numerator / denominator

        if np.isnan(self.estim).any():
            print "Nan in Omega"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Omega"
            raise ValueError


class Beta():
    """
    Class to store and update (M-step) the parameter `beta` in the logistic
    function in the prior of the msCentipede model.

    Arguments
        scores : array
        an array of scores for each motif instance. these could include
        PWM score, conservation score, a measure of various histone
        modifications, outputs from other algorithms, etc.

    """

    def __init__(self, scores):
    
        self.S = scores.shape[1]
        self.estim = np.random.rand(self.S)

    def update(self, scores, zeta):
        """Update the estimates of parameter `beta` in the model.
        """

        def function(x, kwargs):
            """Computes part of the likelihood function that has
            terms containing `beta`.
            """

            scores = kwargs['scores']
            zeta = kwargs['zeta']

            arg = insum(x * scores,[1])
            func = arg * zeta.estim[:,1:] - nplog(1 + np.exp(arg))
            f = -1. * func.sum()
            return f

        def gradient(x, kwargs):
            """Computes gradient of the likelihood function with 
            respect to `beta`.
            """

            scores = kwargs['scores']
            zeta = kwargs['zeta']

            arg = insum(x * scores,[1])
            Df = -1 * np.sum(scores * (zeta.estim[:,1:] - logistic(-arg)),0)
            return Df

        def hessian(x, kwargs):
            """Computes hessian of the likelihood function with 
            respect to `beta`.
            """

            scores = kwargs['scores']
            zeta = kwargs['zeta']

            arg = insum(x * scores,[1])
            larg = scores * logistic(arg) * logistic(-arg)
            Hf = np.dot(scores.T, larg)
            return Hf

        xo = self.estim.copy()
        args = dict([('scores',scores),('zeta',zeta)])
        self.estim = optimizer(xo, function, gradient, hessian, args)

        if np.isnan(self.estim).any():
            print "Nan in Beta"
            raise ValueError

        if np.isinf(self.estim).any():
            print "Inf in Beta"
            raise ValueError        


def optimizer(xo, function, gradient, hessian, kwargs):
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

        # compute likelihood function
        f = function(xx, kwargs)
        if np.isnan(f) or np.isinf(f):
            f = np.array([np.finfo('float32').max]).astype('float')
        else:
            f = np.array([f]).astype('float')
        
        # compute gradient
        Df = gradient(xx, kwargs)
        if np.isnan(Df).any() or np.isinf(Df).any():
            Df = -1 * np.finfo('float32').max * np.ones((1,xx.size), dtype=float)
        else:
            Df = Df.reshape(1,xx.size)
        if z is None:
            return cvx.matrix(f), cvx.matrix(Df)

        # compute hessian
        hess = hessian(xx, kwargs)
        Hf = z[0] * hess
        return cvx.matrix(f), cvx.matrix(Df), cvx.matrix(Hf)

    # warm start for the optimization
    V = xo.size
    x_init = xo.reshape(V,1)

    # call the optimization subroutine in cvxopt
    if kwargs.has_key('G'):
        # call a constrained nonlinear solver
        solution = solvers.cp(F, G=cvx.matrix(kwargs['G']), h=cvx.matrix(kwargs['h']))
    else:
        # call an unconstrained nonlinear solver
        solution = solvers.cp(F)

    x_final = np.array(solution['x']).ravel()

    return x_final


def compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model):
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

    lhood_bound = Data()
    lhood_unbound = Data()

    for j in xrange(data.J):
        value = outsum(data.value[j])[0]
        total = outsum(data.total[j])[0]
        
        lhood_bound.value[j] = outsum([gammaln(data.value[j][r] + pi.value[j] * tau.estim[j]) \
            + gammaln(data.total[j][r] - data.value[j][r] + (1 - pi.value[j]) * tau.estim[j]) \
            - gammaln(data.total[j][r] + tau.estim[j]) + gammaln(tau.estim[j]) \
            - gammaln(pi.value[j] * tau.estim[j]) - gammaln((1 - pi.value[j]) * tau.estim[j]) \
            for r in xrange(data.R)])[0]

        if model in ['msCentipede','msCentipede_flexbgmean']:
            
            lhood_unbound.value[j] = value * nplog(pi_null.value[j]) \
                + (total - value) * nplog(1 - pi_null.value[j])
    
        elif model=='msCentipede_flexbg':
            
            lhood_unbound.value[j] = outsum([gammaln(data.value[j][r] + pi_null.value[j] * tau_null.estim[j]) \
                + gammaln(data.total[j][r] - data.value[j][r] + (1 - pi_null.value[j]) * tau_null.estim[j]) \
                - gammaln(data.total[j][r] + tau_null.estim[j]) + gammaln(tau_null.estim[j]) \
                - gammaln(pi_null.value[j] * tau_null.estim[j]) - gammaln((1 - pi_null.value[j]) * tau_null.estim[j]) \
                for r in xrange(data.R)])[0]

    return lhood_bound, lhood_unbound


def likelihood(data, scores, zeta, pi, tau, \
    alpha, beta, omega, pi_null, tau_null, model):
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

    apriori = insum(beta.estim * scores,[1])

    lhoodA, lhoodB = compute_footprint_likelihood(data, pi, tau, pi_null, tau_null, model)

    footprint = np.zeros((data.N,1),dtype=float)
    for j in xrange(data.J):
        footprint += insum(lhoodA.value[j],[1])

    P_1 = footprint + insum(gammaln(zeta.total + alpha.estim[:,1]) - gammaln(alpha.estim[:,1]) \
        + alpha.estim[:,1] * nplog(omega.estim[:,1]) + zeta.total * nplog(1 - omega.estim[:,1]), [1])
    P_1[P_1==np.inf] = MAX
    P_1[P_1==-np.inf] = -MAX

    null = np.zeros((data.N,1), dtype=float)
    for j in xrange(data.J):
        null += insum(lhoodB.value[j],[1])

    P_0 = null + insum(gammaln(zeta.total + alpha.estim[:,0]) - gammaln(alpha.estim[:,0]) \
        + alpha.estim[:,0] * nplog(omega.estim[:,0]) + zeta.total * nplog(1 - omega.estim[:,0]), [1])
    P_0[P_0==np.inf] = MAX
    P_0[P_0==-np.inf] = -MAX

    L = P_0 * zeta.estim[:,:1] + insum(P_1 * zeta.estim[:,1:],[1]) + apriori * (1 - zeta.estim[:,:1]) \
        - nplog(1 + np.exp(apriori)) - insum(zeta.estim * nplog(zeta.estim),[1])
    
    L = L.sum() / data.N

    if np.isnan(L):
        print "Nan in LogLike"
        return -np.inf

    if np.isinf(L):
        print "Inf in LogLike"
        return -np.inf

    return L


def EM(data, scores, zeta, pi, tau, alpha, beta, omega, pi_null, tau_null, model):
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


    # update binding posteriors
    zeta.update(data, scores, pi, tau, \
            alpha, beta, omega, pi_null, tau_null, model)

    # update multi-scale parameters
    starttime = time.time()
    pi.update(data, zeta, tau)
    print "p_jk update in %.3f secs"%(time.time()-starttime)

    starttime = time.time()
    tau.update(data, zeta, pi)
    print "tau update in %.3f secs"%(time.time()-starttime)

    # update negative binomial parameters
    starttime = time.time()
    omega.update(zeta, alpha)
    print "omega update in %.3f secs"%(time.time()-starttime)

    starttime = time.time()
    alpha.update(zeta, omega)
    print "alpha update in %.3f secs"%(time.time()-starttime)

    # update prior parameters
    starttime = time.time()
    beta.update(scores, zeta)
    print "beta update in %.3f secs"%(time.time()-starttime)

def square_EM(data, scores, zeta, pi, tau, alpha, beta, omega, pi_null, tau_null, model):
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
    # a better likelihood than the best of the two steps.
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


def estimate_optimal_model(reads, totalreads, scores, background, model, restarts, mintol):
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

    # transform data into multiscale representation
    data = Data(reads)
    data_null = Data(background)
    scores = np.hstack((np.ones((data.N,1), dtype=float), scores))
    del reads

    # set background model
    pi_null = Pi(data_null.J)
    for j in xrange(pi_null.J):
        pi_null.value[j] = np.sum(np.sum(data_null.value[j],0),0) / np.sum(np.sum(data_null.total[j],0),0).astype('float')
    tau_null = Tau(data_null.J)
    tau_null = None
        
    if model=='msCentipede_flexbg':
        
        tau_null = Tau(data_null.J)

        zeta_null = Zeta(data_null, background.sum(1))
        zeta_null.estim[:,1] = 1
        zeta_null.estim[:,0] = 0

        # iterative update of background model; 
        # evaluate convergence based on change in estimated
        # background overdispersion
        change = np.inf
        while change>1e-2:
            change = tau_null.estim.copy()
            
            tau_null.update(data_null, zeta_null, pi_null)
            pi_null.update(data_null, zeta_null, tau_null)

            change = np.abs(change-tau_null.estim).sum() / tau_null.J

    maxLoglike = -np.inf
    restart = 0
    err = 1
    runlog = ['Number of sites = %d'%data.N]
    while restart<restarts:

        try:
            totaltime = time.time()
            print "Restart %d ..."%(restart+1)

            # initialize multi-scale model parameters
            pi = Pi(data.J)
            tau = Tau(data.J)

            # initialize negative binomial parameters
            alpha = Alpha(data.R)
            omega = Omega(data.R)

            # initialize prior parameters
            beta = Beta(scores)

            # initialize posterior over latent variables
            zeta = Zeta(data, totalreads)
            for j in xrange(pi.J):
                pi.value[j] = np.sum(data.value[j][0] * zeta.estim[:,1:],0) \
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
            print Loglike

            tol = np.inf
            iter = 0

            while np.abs(tol)>mintol:

                itertime = time.time()
                EM(data, scores, zeta, pi, tau, \
                        alpha, beta, omega, pi_null, tau_null, model)

                newLoglike = likelihood(data, scores, zeta, pi, tau, \
                        alpha, beta, omega, pi_null, tau_null, model)

                tol = newLoglike - Loglike
                Loglike = newLoglike
                print "Iteration %d: log likelihood = %.7f, change in log likelihood = %.7f, iteration time = %.3f secs"%(iter+1, Loglike, tol, time.time()-itertime)
                iter += 1
            totaltime = (time.time()-totaltime)/60.

            # test if mean cleavage rate at bound sites is greater than at 
            # unbound sites, for each replicate; avoids local optima issues.
            negbinmeans = alpha.estim * (1-omega.estim)/omega.estim
            if np.any(negbinmeans[:,0]<negbinmeans[:,1]):
                restart += 1
                log = "%d. Log likelihood (per site) = %.3f (Completed in %.3f minutes)"%(restart,Loglike,totaltime)
                runlog.append(log)
                # choose these parameter estimates, if the likelihood is greater.
                if Loglike>maxLoglike:
                    maxLoglikeres = Loglike
                    if model in ['msCentipede','msCentipede_flexbgmean']:
                        footprint_model = (pi, tau, pi_null)
                    elif model=='msCentipede_flexbg':
                        footprint_model = (pi, tau, pi_null, tau_null)
                    count_model = (alpha, omega)
                    prior = beta

        except ValueError:

            print "encountered an invalid value"
            if err<5:
                print "re-initializing learning for Restart %d ... %d"%(restart,err)
                err += 1
            else:
                print "Error in learning model parameters. Please ensure the inputs are all valid"
                sys.exit(1)

    return footprint_model, count_model, prior, runlog


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

    (N,L,R) = reads.shape
    data = Data(reads)
    data_null = Data(background)
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
        pi_null.value[j] = np.sum(np.sum(data_null.value[j],0),0) \
            / np.sum(np.sum(data_null.total[j],0),0).astype('float')
    tau_null = None

    if model=='msCentipede_flexbg':

        tau_null = footprint[3]

        if data_null.N>1000:

            zeta_null = Zeta(data_null, background.sum(1))
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

    zeta = Zeta(data, totalreads, infer=True)

    zeta.infer(data, scores, pi, tau, alpha, beta, omega, \
        pi_null, tau_null, model)
    
    return zeta.posterior_log_odds, \
        zeta.prior_log_odds, zeta.footprint_log_likelihood_ratio, \
        zeta.total_log_likelihood_ratio

