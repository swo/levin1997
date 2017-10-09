#!/usr/bin/env python3

# Simulation extending Levin 1997. I added a "fraction zero", so that some
# hosts never get treated.

import numpy as np
import sys

class Levin:
    def __init__(self, s, fc, cac, f, g, tg, N, n_generations, H0=1e-9, E0=1e-9, log=None, log_interval=100):
        self.s = s # selection coefficient favoring susceptible strains
        self.fc = fc # fraction consumers
        self.cac = cac # average treatments per year among consumers
        self.f = f # fraction of environment replaced by host-shed bugs
        self.g = g # fraction of each host's bugs replaced by environmental bugs
        self.tg = tg # generation time, in years
        self.N = N # number of hosts
        self.n_generations = n_generations # simulation duration
        self.E = E0 # initial resistance in environment compartment
        self.H = np.full(N, H0) # initial resistance in all hosts
        self.log = log
        self.log_interval = log_interval

        self.n_consumer = int(self.fc * self.N)
        self.n_zero = self.N - self.n_consumer
        self.p_treatment = self.cac * self.tg

        assert 0 <= self.p_treatment <= 1

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
        nonconsumer_treatment = np.full(self.n_zero, 0)
        consumer_treatment = np.random.binomial(1, self.p_treatment, size=self.n_consumer)
        treatment = np.concatenate((nonconsumer_treatment, consumer_treatment))
        self.H = np.fmax(self.H, treatment)

        self.E = self.E + self.f * (np.mean(self.H) - self.E) - self.s * self.E * (1 - self.E) / (1.0 - self.s * self.E)


f = 0.05 # fraction of environment replaced by host-shed bugs
g = 0.005 # fraction of each host's bugs replaced by environmental bugs
tg = 1/219 # generation time, in years
N = 10000 # number of hosts
n_generations = int(5e3)

#ss = [0.04, 0.03, 0.02, 0.01, 0.005]
ss = [0.04, 0.02, 0.005]
cacs = np.linspace(0.0, 2.0, 5)
fnzs = np.linspace(0.0, 1.0, num=5)

print('s', 'fnz', 'cac', 'host', 'env', sep='\t')
for s in ss:
    for cac in cacs:
        for fnz in fnzs:
            host, env = Levin(s, fnz, cac, f, g, tg, N, n_generations).run()
            print(s, fnz, cac, host, env, sep='\t')
            sys.stdout.flush()
