#!/usr/bin/env python3

# Extension to the Levin 1997 simulation. In the original simulation, the resistant strain had
# fitness 1-s without treatment and fitness >0 during treatment. The sensitive strain has fitness
# 1 without treatment and 0 with.
#
# Here, I vary the fitness costs sR in the resistant strain independent, but I also vary the fitness
# 1-sS of the sensitive strain during antibiotic treatment. sS=1 recapitulates the original model.

import numpy as np
import sys

class Levin:
    def __init__(self, sS, sR, fc, cac, f, g, tg, N, n_generations, H0=1e-9, E0=1e-9, log=None, log_interval=100):
        self.sS = sS # fitness cost hurting S strain during abx treatment
        self.sR = sR # fitness cost hurting R strain in environment and untreated hosts
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

    @staticmethod
    def select_fraction(fA, wA, wB):
        '''New fraction fA when strain A has fitness wA and B has wB'''
        return wA * fA / (wA * fA + wB * (1.0 - fA))

    def advance(self):
        # create a vector: 0 means no treatment, 1 means treatment
        nonconsumer_treatment = np.full(self.n_zero, 0)
        consumer_treatment = np.random.binomial(1, self.p_treatment, size=self.n_consumer)
        treatment = np.concatenate((nonconsumer_treatment, consumer_treatment))

        wR = 1.0 - self.sR # strain R has the same fitness in all cases
        wS = 1.0 - self.sS * treatment # strain S has fitness 1 when untreated, but 1-sS when treated

        # prepare to update E using previous H values
        # strain S has fitness 1 in the environment
        new_E = self.f * np.mean(self.H) + (1.0 - self.f) * self.select_fraction(self.E, wR, 1.0)

        # update H with previous E value as if no hosts got drug
        self.H = self.g * self.E + (1.0 - self.g) * self.select_fraction(self.H, wR, wS)

        # update E
        self.E = new_E

f = 0.05 # fraction of environment replaced by host-shed bugs
g = 0.005 # fraction of each host's bugs replaced by environmental bugs
tg = 1/219 # generation time, in years
N = 10000 # number of hosts
n_generations = int(5e3)

sSs = np.linspace(0.75, 1.0, num=4)
sRs = [0.04, 0.02, 0.005]
cacs = np.linspace(0.0, 2.0, num=5)
fnzs = np.linspace(0.0, 1.0, num=5)

print('sS', 'sR', 'fnz', 'cac', 'host', 'env', sep='\t')
for sS in sSs:
    for sR in sRs:
        for cac in cacs:
            for fnz in fnzs:
                host, env = Levin(sS, sR, fnz, cac, f, g, tg, N, n_generations).run()
                print(sS, sR, fnz, cac, host, env, sep='\t')
                sys.stdout.flush()
