#!/usr/bin/env python
import sys
import operator
from math import log, lgamma
from collections import Counter


######################################### Likelihood part #########################################

def cluster(cluster1, cluster2):
    return cluster_likelihood(len(cluster1.subjects) + len(cluster2.subjects),
                              cluster1.feature_counts + cluster2.feature_counts)

def cluster_likelihood(num_subjects, feature_counter):
    if len(feature_sizes) == 0 and len(feature_prior_counts) == 0:
        feature_count_freqs = Counter(feature_counter.values())
        n0 = num_features - len(feature_counter)
        likelihood = n0 * cluster_likel(num_subjects, 0) if n0 > 0 else 0
        for x in feature_count_freqs:
            likelihood += feature_count_freqs[x] * cluster_likel(num_subjects, x)
        return likelihood

    likelihood = 0
    n0 = num_features
    for f in feature_counter:
        feature_cluster_size = 1 if f not in feature_sizes else feature_sizes[f]
        feature_cluster_prior_count = int(num_subjects * feature_cluster_size / 2.0) if f not in feature_prior_counts else feature_prior_counts[f]
        p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_cluster_size)
        feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
        likelihood += cluster_likel(num_subjects * feature_cluster_size, feature_counter[f], feature_cluster_alpha, feature_cluster_beta)
        n0 -= feature_cluster_size
    for f in feature_sizes:
        if f not in feature_counter:
            feature_cluster_prior_count = int(num_subjects / 2.0) if f not in feature_prior_counts else feature_prior_counts[f]
            p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_sizes[f])
            feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)

            likelihood += cluster_likel(num_subjects * feature_sizes[f], 0, feature_cluster_alpha, feature_cluster_beta)
            n0 -= feature_sizes[f]
        else:
            for f in feature_prior_counts:
                if f not in feature_counter and f not in feature_sizes:
                    p = (feature_prior_counts[f] + 1) / float((total_num_subjects + 2))
                    feature_cluster_alpha,feature_cluster_beta = beta_parameters(num_features, p)
                    likelihood += cluster_likel(num_subjects * feature_prior_counts[f], 0, feature_cluster_alpha, feature_cluster_beta)
                    n0 - 1
    if n0 > 0:
        feature_cluster_prior_count = int(num_subjects * feature_cluster_size / 2.0)
        p = (feature_cluster_prior_count + 1) / float((total_num_subjects + 2) * feature_cluster_size)
        feature_cluster_alpha, feature_cluster_beta = beta_parameters(num_features, p)
        likelihood += n0 * cluster_likel(num_subjects, 0, feature_cluster_alpha, feature_cluster_beta)
            
    return likelihood


def beta_parameters(n, mu):
     s = ((mu * (1 - mu)) / float(n))**0.5
     a = ((1 - mu) / float(s)- 1 / float(mu)) * (mu**2)
     b = a * (1 / float(mu) - 1)
     return a, b


def cluster_likel(n, x, alpha=1, beta=1):
  if x > n:
      sys.exit("Error: x>n\n")
  return lgamma(alpha + beta) + lgamma(alpha + x) + lgamma(beta + n - x) - lgamma(alpha) -\
         lgamma(beta) - lgamma(alpha + beta + n)



########################################################################################
###############################   CLASS CLUSTER     ####################################
########################################################################################

class Cluster:
  count = 0
  lik = 0
  def __init__(self, subject, features=[]):
    try:
        self.subjects = [int(subject)]
    except:
        self.subjects = [subject]
    if all([":" in f for f in features]):
        self.feature_counts = Counter()
        for f in features:
            key_val = f.split(":")
            self.feature_counts[key_val[0]] += int(key_val[1])
    else:
        self.feature_counts = Counter(features)
    self.cluster_similarity = {}
    Cluster.count += 1

  def close(self):
    del self.subjects
    del self.feature_counts
    del self.lik

  def set_lik(self):
    self.lik = cluster_likelihood(len(self.subjects), self.feature_counts)

  def merge(self, other_cluster, sim=""):
    self.subjects += other_cluster.subjects
    self.feature_counts += other_cluster.feature_counts
    self.lik = sim + self.lik + other_cluster.lik
    other_cluster.close()
    del other_cluster
    Cluster.count = Cluster.count - 1

####################################################################### end of class

def initialise():
 global total_num_subjects, num_features, feature_sizes, feature_prior_counts, informative_prior
 feature_separator = ", "
 cluster_list = []
 feature_totals = Counter()
 feature_sizes = {}
 feature_prior_counts = {}
 informative_prior = True
 try:
     F = open(sys.argv[1], "r")
     for line in F:
         d = line.strip().split()
         feature_sizes[d[0]] = int(d[1])
 except:
     None
 total_num_subjects = 0
 for line in sys.stdin:
     key_val = line.strip().replace("U","").replace("C","").split("\t")
     try:
         cluster_list += [Cluster(key_val[0], key_val[1].split(feature_separator))]
         try:
             total_num_subjects = max(total_num_subjects, int(key_val[0]))

         except:
             None
     except:
         cluster_list += [Cluster(Cluster.count,key_val[0].split(feature_separator))]
     for feature in key_val[-1].split(feature_separator):
         feature_totals[feature] += 1
 total_num_subjects = max(total_num_subjects, len(cluster_list))
 try:
     F = open(sys.argv[2], "r")
     for line in F:
         d = line.strip().split()
         feature_prior_counts[d[0]] = int(d[1])
     informative_prior = True
 except:
     if informative_prior:
         feature_prior_counts = feature_totals
 num_features = len(feature_totals)
 Cluster.likel = 0
 for cl in cluster_list:
     cl.set_lik()
     Cluster.lik += cl.lik
 for i in range(Cluster.count - 1):
     for j in range(i + 1, Cluster.count):
         sim = cluster_similarity(cluster_list[i], cluster_list[j])
         cluster_list[i].cluster_similarity[cluster_list[j]] = sim
         cluster_list[j].cluster_similarity[cluster_list[i]] = sim
 return cluster_list


def cluster_similarity(cluster1, cluster2):
    sim = cluster(cluster1, cluster2) - cluster1.lik - cluster2.lik
    return sim


def find_most_similar():
  max_similarity= -sys.float_info.max
  for i in range(Cluster.count - 1):
    for j in range(i + 1, Cluster.count):
      similarity = cluster_list[i].cluster_similarity[cluster_list[j]]
      if similarity > max_similarity:
        best_pair = [i,j]
        max_similarity = similarity
  return best_pair, max_similarity


def merge_clusters(i, j, sim=""):
  if i > j:
      temp = i
      i = j
      j = temp
  for cl in cluster_list:
      if cl != cluster_list[i] and cl != cluster_list[j]:
          cl.cluster_similarity.pop(cluster_list[j])
  cluster_list[i].merge(cluster_list[j], sim)
  cluster_list.pop(j)
  if i > j:
      i = i - 1
  for cl in cluster_list:
      if cl != cluster_list[i]:
          sim = cluster_similarity(cluster_list[i], cl)
          cluster_list[i].cluster_similarity[cl] = sim
          cl.cluster_similarity[cluster_list[i]] = sim


def get_clustering_from_cluster_list():
    clustering = {}
    counter = 0
    for cl in cluster_list:
        counter = counter + 1
        for sub in cl.subjects:
            clustering[sub] = counter
    return clustering


def user_clustering():
  max_l = Cluster.lik
  while Cluster.count > 1:
      pair, delta = find_most_similar()
      if delta > 0:
          merge_clusters(pair[0], pair[1], delta)
          Cluster.lik += delta
          if Cluster.lik > max_l:
              max_l = Cluster.lik
      else:
          break
  return get_clustering_from_cluster_list()


##################################################################################
######################################## MAIN ####################################

def main():
  global cluster_list, c, y
  cluster_list = initialise()
  clustering = user_clustering()
  print(sorted(clustering.items(), key=operator.itemgetter(0)))


if __name__ == '__main__':
    main()





    
