#!/usr/bin/env python3

# Simulation described in Levin, et al. "Population genetics of antibiotic resistance",
# CID 1997.

from levin import sim

s = 0.01 # selection coefficient favoring susceptible strains
T = 0.8 # average treatments per year
f = 0.05 # fraction of environment replaced by host-shed bugs
g = 0.005 # fraction of each host's bugs replaced by environmental bugs
tg = 1/219 # generation time, in years
N = 10000 # number of hosts

print(sim(s, T, f, g, tg, N, int(5e2)))
