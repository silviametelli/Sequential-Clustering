#!/usr/bin/env python
import sys
import random
import scipy
import scipy.stats
from scipy.stats import bernoulli,uniform,poisson

n=k_n=100 if len(sys.argv)<2 else int(sys.argv[1])
p=k_p=100 if len(sys.argv)<3 else int(sys.argv[2])
mean_n_cluster_size=10 if len(sys.argv)<4 else int(sys.argv[3])
mean_p_cluster_size=10 if len(sys.argv)<5 else int(sys.argv[4])
while k_n>=n:
    k_n=poisson.rvs(n/float(mean_n_cluster_size))
while k_p>=p:
    k_p=poisson.rvs(p/float(mean_p_cluster_size))
n_breaks=[0]+sorted(random.sample(range(1,n),k_n))+[n+1]
p_breaks=[0]+sorted(random.sample(range(1,p),k_p))+[p+1]
bernoulli_pars=uniform.rvs(size=(k_p+1)*(k_n+1))
for c_n in range(k_n+1):
    for c_p in range(k_p+1):
        theta=bernoulli_pars[c_n*(k_p+1)+c_p]
        for i in range(n_breaks[c_n],n_breaks[c_n+1]):
            for j in range(p_breaks[c_p],p_breaks[c_p+1]):
                if bernoulli.rvs(theta):
                    print str(i+1)+","+str(j+1)
