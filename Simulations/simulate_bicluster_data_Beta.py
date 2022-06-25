#!/usr/bin/env python
import sys
import random
from scipy.stats import bernoulli, binom, poisson, beta

p_cluster_counts = False
default_kn, default_kp = 100, 100
n = k_n = default_kn if len(sys.argv) < 2 else int(sys.argv[1])
p = k_p = default_kp if len(sys.argv) < 3 else int(sys.argv[2])
mean_n_cluster_size = 10 if len(sys.argv) < 4 else int(sys.argv[3])
mean_p_cluster_size = 10 if len(sys.argv) < 5 else int(sys.argv[4])
while k_n >= n:
    k_n = poisson.rvs(n / float(mean_n_cluster_size))
while k_p >= p:
    k_p=poisson.rvs(p / float(mean_p_cluster_size))
n_breaks = [0] + sorted(random.sample(range(1, n), k_n)) + [n]
p_breaks = [0] + sorted(random.sample(range(1, p), k_p)) + [p]
bernoulli_pars = [beta.rvs(0.1, 1) if bernoulli.rvs(.5) else beta.rvs(1,0.1) for _ in range((k_p + 1) * (k_n + 1))]
for c_n in range(k_n + 1):
    if not p_cluster_counts:
        for c_p in range(k_p + 1):
            theta = bernoulli_pars[c_n * (k_p + 1) + c_p]
            for i in range(n_breaks[c_n], n_breaks[c_n + 1]):
                for j in range(p_breaks[c_p], p_breaks[c_p + 1]):
                    if bernoulli.rvs(theta):
                        print(str(i + 1) + "," + str(j + 1))
    else:
        for i in range(n_breaks[c_n], n_breaks[c_n + 1]):
            n_i_counts = [None] * (k_p + 1)
            for c_p in range(k_p + 1):
                theta = bernoulli_pars[c_n * (k_p + 1) + c_p]
                n_i_counts[c_p] = binom.rvs(p_breaks[c_p + 1] - p_breaks[c_p], theta)
            non_zero_elements = [c_p for c_p in range(k_p + 1) if n_i_counts[c_p] > 0]

            if len(non_zero_elements) > 0:
                print (str(i + 1) + "\t" + ", ".join([str(c_p + 1)+ ":" + str(n_i_counts[c_p]) for c_p in non_zero_elements]))
