#!/usr/bin/env python3

# Simulation described in Levin, et al. "Population genetics of antibiotic resistance",
# CID 1997.

import numpy as np
import sys

class Levin:
    def __init__(self, s, T, f, g, tg, N, n_generations, H0=1e-9, E0=1e-9, log=None, log_interval=100):
        self.s = s
        self.T = T
        self.f = f
        self.g = g
        self.tg = tg
        self.N = N
        self.n_generations = n_generations
        self.E = E0
        self.H = np.full(N, H0)
        self.log = log
        self.log_interval = log_interval

        self.p_treatment = T * tg

    def run(self):
        if self.log is not None:
            # print the header
            print('generation', 'host', 'env', sep='\t', file=self.log)

            n_complete = 0
            while n_complete < self.n_generations:
                this_interval = min(self.log_interval, self.n_generations - n_complete)
                self.advance_n(this_interval)
                n_complete += this_interval
                print(n_complete, np.mean(self.H), self.E, sep='\t', file=self.log)
        else:
            self.advance_n(self.n_generations)

        return (np.mean(self.H), self.E)

    def advance_n(self, n_generations):
        for i in range(n_generations):
            self.advance()

    def advance(self):
        # act as if no hosts got drug
        self.H = self.H + self.g * (self.E - self.H) - self.s * self.H * (1.0 - self.H) / (1.0 - self.s * self.H)

        # then replace some with treatment
        treatment = np.random.binomial(1, self.p_treatment, size=self.N).astype(float)
        self.H = np.fmax(self.H, treatment)

        self.E = self.E + self.f * (np.mean(self.H) - self.E) - self.s * self.E * (1 - self.E) / (1.0 - self.s * self.E)


s = 0.005 # selection coefficient favoring susceptible strains
T = 2.0 # average treatments per year
f = 0.05 # fraction of environment replaced by host-shed bugs
g = 0.005 # fraction of each host's bugs replaced by environmental bugs
tg = 1/219 # generation time, in years
N = 10000 # number of hosts
n_generations = int(5e3)

print('s', 'T', 'host', 'env', sep='\t')
for s in [0.04, 0.03, 0.02, 0.01, 0.005]:
    for T in np.arange(0, 2.01, 0.2):
        host, env = Levin(s, T, f, g, tg, N, n_generations).run()
        print(s, T, host, env, sep='\t')
        sys.stdout.flush()
