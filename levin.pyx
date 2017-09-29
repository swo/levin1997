# Simulation described in Levin, et al. "Population genetics of antibiotic resistance",
# CID 1997.

import numpy as np
cimport numpy as np
DTYPE = np.float64

def sim(float s, T, float f, g, tg, int N, int n_generations, E0=1e-9, H0=1e-9):
    '''
    s : float
      selection coefficient favoring susceptible strains
    T : float
      average treatments per year
    f : float
      fraction of environment replaced by host-shed bugs
    g : float
      fraction of each host's bugs replaced by environmental bugs
    tg : float
      generation time, in years
    N : int
      number of hosts
    n_generations : int
      number of generations
    E0 : float
      initial fraction resistant in environment
    H0 : float
      initial fraction resistant in hosts
    '''

    cdef int generation_i, host_i
    cdef float p_treatment, E

    p_treatment = T * tg # p_treatment * 1 year / (tx per year) = T

    E = E0
    H = np.full(N, 1e-9, DTYPE)

    for generation_i in range(n_generations):
        # compute new environment
        # environment gets host bugs and experiences in-environment selection
        E_new = E + f * (np.mean(H) - E) - s * E * (1 - E) / (1.0 - s * E)

        # compute & update hosts
        for host_i in range(N):
            # does this host get treatment?
            if np.random.binomial(1, p_treatment):
                # treated hosts get 100% resistance
                H[host_i] = 1.0
            else:
                # untreated hosts get environmental bugs and experience in-host selection
                H[host_i] = H[host_i] + g * (E - H[host_i]) - s * H[host_i] * (1.0 - H[host_i]) / (1.0 - s * H[host_i])

        # update environment
        E = E_new

    return (np.mean(H), E)
