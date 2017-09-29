#!/usr/bin/env python3

# Simulation described in Levin, et al. "Population genetics of antibiotic resistance",
# CID 1997.

from numba import jit
import numpy as np

@jit
def sim(s, T, f, g, tg, N, n_generations):
    p_treatment = T * tg # p_treatment * 1 year / (tx per year) = T

    E = 1e-9 # initial fraction resistant in environment
    H = np.array([1e-9 for i in range(N)]) # initial fraction resistant in each host

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

    # make report
    return np.mean(H)

s = 0.01 # selection coefficient favoring susceptible strains
T = 2.0 # average treatments per year
f = 0.05 # fraction of environment replaced by host-shed bugs
g = 0.005 # fraction of each host's bugs replaced by environmental bugs
tg = 1/219 # generation time, in years
N = 10000 # number of hosts
n_generations = int(1e4)

sim(s, T, f, g, tg, N, n_generations)
