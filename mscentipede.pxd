import numpy as np
cimport numpy as np
from cpython cimport bool

cdef class Data:

	cdef public long N, L, R, J
	cdef public dict valueA, valueB, total

	cdef transform_to_multiscale(self, np.ndarray[np.float64_t, ndim=3] reads)


cdef class Zeta:

	cdef public long N
	cdef public np.ndarray total, prior_log_odds, \
		footprint_log_likelihood_ratio, total_log_likelihood_ratio, \
		posterior_log_odds, estim

	cdef update(self, Data data, np.ndarray[np.float64_t, ndim=2] scores, \
        Pi pi, Tau tau, Alpha alpha, Beta beta, Omega omega, \
        Pi pi_null, Tau tau_null, str model)

	cdef infer(self, Data data, np.ndarray[np.float64_t, ndim=2] scores, \
        Pi pi, Tau tau, Alpha alpha, Beta beta, Omega omega, \
        Pi pi_null, Tau tau_null, str model)


cdef class Pi:

	cdef public long J
	cdef public dict value

cpdef tuple pi_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args)

cpdef tuple pi_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args)

cdef class Tau:

	cdef public long J
	cdef public np.ndarray estim

cpdef tuple tau_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args)

cpdef tuple tau_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args)

cdef class Alpha:

	cdef public long R
	cdef public np.ndarray estim

cpdef tuple alpha_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args)

cpdef tuple alpha_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args)

cdef class Omega:

	cdef public long R
	cdef public np.ndarray estim

	cdef update(self, Zeta zeta, Alpha alpha)


cdef class Beta:

	cdef public long S
	cdef public np.ndarray estim

cpdef tuple beta_function_gradient(np.ndarray[np.float64_t, ndim=1] x, dict args)

cpdef tuple beta_function_gradient_hessian(np.ndarray[np.float64_t, ndim=1] x, dict args)

cdef tuple compute_footprint_likelihood(Data data, Pi pi, Tau tau, Pi pi_null, Tau tau_null, str model)


cdef double likelihood(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
	Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
	Omega omega, Pi pi_null, Tau tau_null, str model)

cdef EM(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
    Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
    Omega omega, Pi pi_null, Tau tau_null, str model)

cdef square_EM(Data data, np.ndarray[np.float64_t, ndim=2] scores, \
    Zeta zeta, Pi pi, Tau tau, Alpha alpha, Beta beta, \
    Omega omega, Pi pi_null, Tau tau_null, str model)
