#!/usr/bin/env python
import sys
import operator
from math import log, lgamma
import timeit
from collections import Counter



#start = timeit.default_timer()

def cluster(cluster1,cluster2):
    return cluster_likelihood(len(cluster1.subjects)+len(cluster2.subjects),cluster1.feature_counts+cluster2.feature_counts)
 
def cluster_likelihood(num_subjects,feature_counter):
    if len(feature_sizes)==0:
        feature_count_freqs=Counter(feature_counter.values())
        n0=num_features-len(feature_counter)
        likelihood=n0*cluster_likel(num_subjects,0) if n0>0 else 0
        for x in feature_count_freqs:
            likelihood+=feature_count_freqs[x]*cluster_likel(num_subjects,x)
        return likelihood

    likelihood=0
    n0=num_features
    for f in feature_counter:
        feature_cluster_size=1 if f not in feature_sizes else feature_sizes[f]
        likelihood+=cluster_likel(num_subjects*feature_cluster_size,feature_counter[f])
        n0-=feature_cluster_size
    for f in feature_sizes:
        if f not in feature_counter:
            likelihood+=cluster_likel(num_subjects*feature_sizes[f],0)
            n0-=feature_sizes[f]
    if n0>0:
        likelihood+=n0*cluster_likel(num_subjects,0)
    return likelihood

def cluster_likel(n,x):
  if x>n:
      sys.stderr.write("Error: x>n\n")
      exit()
  return lgamma(alpha+beta)+lgamma(alpha+x)+lgamma(beta+n-x)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+n)

class Cluster:
  count=0
  lik=0      
  def __init__(self, subject, features=[]):
    try:
        self.subjects=[int(subject)]
    except:
        self.subjects=[subject]
    if all([":" in f for f in features]):
        self.feature_counts=Counter()
        for f in features:
            key_val=f.split(":")
            self.feature_counts[key_val[0]]+=int(key_val[1])
    else:
        self.feature_counts=Counter(features)
    self.cluster_similarity={}
    Cluster.count+=1
  def close(self):
    del self.subjects
    del self.feature_counts
    del self.lik  
  def set_lik(self):
    self.lik=cluster_likelihood(len(self.subjects),self.feature_counts)
  def merge(self, other_cluster,sim=""):
    self.subjects+=other_cluster.subjects
#      self.subjects.sort()
    self.feature_counts+=other_cluster.feature_counts
    self.lik=sim+self.lik+other_cluster.lik
    other_cluster.close()
    del other_cluster
    Cluster.count=Cluster.count-1

def initialise():
 global num_subjects,num_features,feature_sizes
 feature_separator=", "
 cluster_list=[]
 feature_set=set()
 cl_count=0
 feature_sizes={}
 try:
     F=open(sys.argv[1],"r")
     for line in F:
         d=line.strip().split()
         feature_sizes[d[0]]=int(d[1])
 except:
     None
 num_subjects=0
 for line in sys.stdin:
     cl_count+=1
     key_val=line.strip().translate(None,"UC").split("\t")
     try:
         cluster_list+=[Cluster(key_val[0],key_val[1].split(feature_separator))]
         try:
             num_subjects=max(num_subjects,int(key_val[0]))
         except:
             None
     except:
         cluster_list+=[Cluster(cl_count,key_val[0].split(feature_separator))]
     for feature in key_val[-1].split():
         feature_set.add(feature)

 num_subjects=max(num_subjects,len(cluster_list))
 num_features=len(feature_set)
 Cluster.likel=0
 for cl in cluster_list:
     cl.set_lik()
     Cluster.lik+=cl.lik
 for i in range(Cluster.count-1):
     for j in range(i+1, Cluster.count):
         sim=cluster_similarity(cluster_list[i],cluster_list[j])
         cluster_list[i].cluster_similarity[cluster_list[j]]=sim
         cluster_list[j].cluster_similarity[cluster_list[i]]=sim
 return cluster_list

def cluster_similarity(cluster1, cluster2):
    sim=cluster(cluster1,cluster2)-cluster1.lik-cluster2.lik
    return sim

def find_most_similar():
  max_similarity= -sys.float_info.max
  for i in range(Cluster.count-1):
    for j in range(i+1, Cluster.count):
      similarity=cluster_list[i].cluster_similarity[cluster_list[j]]#similarity(cluster_list[i],cluster_list[j])
      if similarity>max_similarity:
        best_pair=[i,j]
        max_similarity=similarity
  return best_pair, max_similarity

def merge_clusters(i,j, sim=""):
  if i>j:
      temp=i
      i=j
      j=temp
  for cl in cluster_list:
      if cl!=cluster_list[i] and cl!=cluster_list[j]:
          cl.cluster_similarity.pop(cluster_list[j])
  cluster_list[i].merge(cluster_list[j], sim)
  cluster_list.pop(j)
  if i>j:
      i=i-1
  for cl in cluster_list:
      if cl!= cluster_list[i]:
          sim=cluster_similarity(cluster_list[i],cl)
          cluster_list[i].cluster_similarity[cl]=sim
          cl.cluster_similarity[cluster_list[i]]=sim
          
def get_clustering_from_cluster_list():
    clustering={}
    counter=0
    for cl in cluster_list:
        counter=counter+1
        for sub in cl.subjects:
            clustering[sub]=counter
    return clustering

def user_clustering():
  max_l=Cluster.lik
  while Cluster.count>1:  
      pair, delta=find_most_similar()
      if delta>0:
 #         print cluster_list[pair[0]].subjects,cluster_list[pair[1]].subjects
          merge_clusters(pair[0],pair[1],delta)
          Cluster.lik+=delta
          if Cluster.lik>max_l:
              max_l=Cluster.lik
          if Cluster.count==1:
              continue_merging=False
      else:
          break
  return get_clustering_from_cluster_list()

def main():
  global data,cluster_list,c,y,alpha,beta
  alpha=1
  beta=0.5
#  data=Data()
  cluster_list=initialise()
#  c=counting_initialise()
  clustering=user_clustering()
  #print "".join(str(l+1) for l in clustering)
  print sorted(clustering.items(),key=operator.itemgetter(0))

#  stop = timeit.default_timer()
# print stop - start 

  

main()





    
