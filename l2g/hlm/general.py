import numpy as np

def single_exponential_qpar(drsep: np.ndarray, Bpm: float, Rb: float,
    lambda_q : float, q_parallel: float) -> np.ndarray:
    r"""Applies a parallel single exponential profile to the provided drsep.

    For the amplitude of the exponential function a simple constant is taken
    and applied.

    .. math::

       q = q_{\parallel} e^{\frac{-drsep}{\lambda_q}}

    Returns:
        q (arr): 1D array of parallel heat loads.

    Arguments:
        drsep
    """
    q = q_parallel * np.exp(-drsep / lambda_q)
    return q

def single_exponential_psol(drsep: np.ndarray, Bt: float, Bpm: float, Rb: float,
    lambda_q: float, P_sol: float, F: float = 0.5) -> np.ndarray:
    r"""Applies a parallel single exponential profile to the provided drsep.

    The amplitude of the exponential function is calculated via the following
    question:

    .. math::

       q_{\parallel} = \frac{F P_{sol} B_t}{2 \pi B_{pm} R_{b}}

    .. math::

       q = q_{\parallel} e^{\frac{-drsep}{\lambda_q}}

    .. note::

       Note that in the equation we are missing the B toroidal at midplane,
       coming as pitch or ratio between poloidal and toroidal. The reason is,
       since if we apply the total flux expansion the B toroidal at midplane
       cancels out.

    Returns:
        q (arr): 1D array of parallel heat loads.

    Arguments:
        drsep (arr): 1D array of distances from the midplane. In meters.

    """

    q_parallel = F * P_sol * Bt / (2 * np.pi * Bpm * Rb * lambda_q)
    q = q_parallel * np.exp(-drsep / lambda_q)

    return q

def double_exponential_psol(drsep: np.ndarray, Bt: float, Bpm: float, Rb: float,
    lambda_q_main: float, lambda_q_near: float, Rq: float, P_sol: float,
    F: float = 0.5) -> np.ndarray:
    r"""The double exponential plasma profile. Not to be confused of the sum
    of two single exponentials.

    .. math:

        q_{\parallel} = \frac{F P_{sol} B_t}{2 \pi R_b B_{pm} (\lambda_{q, main} + R_q \lambda_{q, near})}

        q = q_{\parallel} (e^{\frac{-drsep}{\lambda{q, main}}} + R_q e^{\frac{-drsep}{\lambda{q, near}}})

    """

    q_parallel = F * P_sol * Bt / (2 * np.pi * Bpm * Rb * (lambda_q_main + Rq * lambda_q_near))

    q = q_parallel * (np.exp(-drsep / lambda_q_main) + Rq * np.exp(-drsep / lambda_q_near))

    return q

def double_exponential_qpar(*args, **kwargs):
    """Not implemented yet."""
    raise NotImplementedError
