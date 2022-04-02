#!/usr/bin/env python
import sys
import operator
from math import log, lgamma
from collections import Counter
import linecache




def New_clustering(info,M,Reading_file):
    global last_line
    global cluster_list,feature_totals
    last_line=M
    print 'SEQUENTIAL CLUSTERING: considering line ',M+1


    cluster_list_old=info[0]
    feature_totals=info[1]


    L=len(cluster_list_old)
    cluster_list=initialise2(cluster_list_old,Reading_file)


    # Find highest similarities of all clusters
    # with the last one added
    max_similarity= -sys.float_info.max
    for i in range(L):
        sims_with_last=cluster_list[i].cluster_similarity[cluster_list[-1]]
        # print 'ith&last',cluster(cluster_list[i],cluster_list[-1])
        # print 'self-',i,cluster_list[i].lik
        # print 'self-last',cluster_list[-1].lik

        if sims_with_last>max_similarity:
            max_similarity=sims_with_last
            best_pair=[i,L]
#merging of clusters
    if max_similarity>0:
     cluster_list[best_pair[0]].merge(cluster_list[best_pair[1]], max_similarity)
     cluster_list.pop(L)
     # for i in range(0,len(cluster_list)):
     #    print 'Cluster %d:' %(i+1),cluster_list[i]

    info[0]=cluster_list
    info[1]=feature_totals
    return info,M+1







###############################   INITIALISE     ####################################

def initialise2(cluster_list_old,Reading_file):
 global total_num_subjects,num_features,feature_sizes
 global feature_prior_counts,informative_prior,feature_totals
 global last_line
 feature_separator=","
 cluster_list=cluster_list_old
 # feature_totals=Counter()
 feature_sizes={}
 feature_prior_counts={}
 informative_prior=False
 total_num_subjects=last_line+1
 #read last_line+1
 line=linecache.getline(Reading_file, last_line+1)
 #create cluster with data from line las_line+1
 key_val=line.strip().translate(None,"UC").split("\t")

 try:
     cluster_list+=[Cluster(key_val[0],key_val[1].split(feature_separator))]
   #  try:
    #     total_num_subjects=max(total_num_subjects,int(key_val[0]))
    # except:
         #None
 except:
     cluster_list+=[Cluster(Cluster.count,key_val[0].split(feature_separator))]
 for feature in key_val[-1].split(feature_separator):
      feature_totals[feature]+=1


 #total_num_subjects=max(total_num_subjects,len(cluster_list))

 # if informative_prior:
 #         feature_prior_counts=feature_totals
 #         print feature_prior_counts
 num_features=len(feature_totals)
 # print 'num_features',num_features
 Cluster.lik=0
 for cl in cluster_list:
 #?????? WHY DOES IT NOT WORK WITH THE CLASS FUNCTION SET_LIK()?????????
     # cl.set_lik()
     selfLik=cluster_likelihood(len(cl.subjects),cl.feature_counts)
     cl.lik=selfLik
     Cluster.lik+=cl.lik





 for i in range(len(cluster_list)-1):
         sim=cluster_similarity(cluster_list[i],cluster_list[-1])
         cluster_list[i].cluster_similarity[cluster_list[-1]]=sim
         cluster_list[-1].cluster_similarity[cluster_list[i]]=sim
         # print 'sim',i,'last=',sim
 return cluster_list






###############################   Class CLUSTER     ####################################


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
    print '------We are inside set_lik(self):'
    self.lik=cluster_likelihood(len(self.subjects),self.feature_counts)

  def merge(self, other_cluster,sim=""):
 #   print "Hi, it's merge()"
#    for i in range(0,2):
 #       print cluster_list[i].subjects,
#        print ''
    self.subjects+=other_cluster.subjects
#      self.subjects.sort()
    self.feature_counts+=other_cluster.feature_counts
    self.lik=sim+self.lik+other_cluster.lik
    other_cluster.close()
    del other_cluster
    Cluster.count=Cluster.count-1

  def merge_in(self, other_cluster):

 #     print "Hi, it's merge_in()"
      self.subjects+=other_cluster.subjects
#      self.subjects.sort()
      self.feature_counts+=other_cluster.feature_counts
      self.lik=self.lik+other_cluster.lik
      other_cluster.close()
      del other_cluster
      Cluster.count=Cluster.count-1

  def __str__(self):
    print_string='Cluster[Subjects:'
    for i in range(0,len(self.subjects)):
        print_string += '%d, ' % self.subjects[i]
    print_string += 'Features: '
    print_string += str(list(self.feature_counts.items()))
    print_string+=']'
    return print_string







###############################   CLUSTER(CLUSTER1,CLUSTER2)     ####################################


def cluster(cluster1,cluster2):
    return cluster_likelihood(len(cluster1.subjects)+len(cluster2.subjects),cluster1.feature_counts+cluster2.feature_counts)



def cluster_likelihood(num_subjects,feature_counter):
    if len(feature_sizes)==0 and len(feature_prior_counts)==0:
        feature_count_freqs=Counter(feature_counter.values())
        # print 'feature_count_freqs', feature_count_freqs
        n0=num_features-len(feature_counter)
        likelihood=n0*cluster_likel(num_subjects,0) if n0>0 else 0
        for x in feature_count_freqs:
            likelihood+=feature_count_freqs[x]*cluster_likel(num_subjects,x)
        return likelihood

    likelihood=0
    n0=num_features
    for f in feature_counter:
        feature_cluster_size=1 if f not in feature_sizes else feature_sizes[f]
        feature_cluster_prior_count=int(num_subjects*feature_cluster_size/2.0) if f not in feature_prior_counts else feature_prior_counts[f]
        p=(feature_cluster_prior_count+1)/float((total_num_subjects+2)*feature_cluster_size)
        feature_cluster_alpha,feature_cluster_beta=beta_parameters(num_features,p)
#        feature_cluster_alpha=(feature_cluster_prior_count+1)/float((num_subjects+2)*feature_cluster_size)
#        feature_cluster_beta=(num_subjects*feature_cluster_size-feature_cluster_prior_count+1)/float((num_subjects+2)*feature_cluster_size)
        likelihood+=cluster_likel(num_subjects*feature_cluster_size,feature_counter[f],feature_cluster_alpha,feature_cluster_beta)
        n0-=feature_cluster_size
    for f in feature_sizes:
        if f not in feature_counter:
            feature_cluster_prior_count=int(num_subjects/2.0) if f not in feature_prior_counts else feature_prior_counts[f]
#            feature_cluster_alpha=(feature_cluster_prior_count+1)/float((num_subjects+2)*feature_sizes[f])
#            feature_cluster_beta=(num_subjects*feature_sizes[f]-feature_cluster_prior_count+1)/float((num_subjects+2)*feature_sizes[f])
            p=(feature_cluster_prior_count+1)/float((total_num_subjects+2)*feature_sizes[f])
            feature_cluster_alpha,feature_cluster_beta=beta_parameters(num_features,p)

            likelihood+=cluster_likel(num_subjects*feature_sizes[f],0, feature_cluster_alpha,feature_cluster_beta)
            n0-=feature_sizes[f]
        else:
            for f in feature_prior_counts:
                if f not in feature_counter and f not in feature_sizes:
#                    feature_cluster_alpha=(feature_prior_counts[f]+1)/float((num_subjects+2))
#                    feature_cluster_beta=(num_subjects-feature_prior_counts[f]+1)/float((num_subjects+2))
                    p=(feature_prior_counts[f]+1)/float((total_num_subjects+2))
                    feature_cluster_alpha,feature_cluster_beta=beta_parameters(num_features,p)
                    likelihood+=cluster_likel(num_subjects*feature_prior_counts[f],0,feature_cluster_alpha,feature_cluster_beta)
                    n0-1
    if n0>0:
        feature_cluster_prior_count=int(num_subjects*feature_cluster_size/2.0)
 #       feature_cluster_alpha=(feature_cluster_prior_count+1)/float((num_subjects+2)*feature_cluster_size)
#        feature_cluster_beta=(num_subjects*feature_cluster_size-feature_cluster_prior_count+1)/float((num_subjects+2)*feature_cluster_size)
        p=(feature_cluster_prior_count+1)/float((total_num_subjects+2)*feature_cluster_size)
        feature_cluster_alpha,feature_cluster_beta=beta_parameters(num_features,p)
        likelihood+=n0*cluster_likel(num_subjects,0,feature_cluster_alpha,feature_cluster_beta)
    
    return likelihood







def cluster_likel(n,x,alpha=1,beta=1):
  if x>n:
      sys.stderr.write("Error: x>n\n")
      exit()
  return lgamma(alpha+beta)+lgamma(alpha+x)+lgamma(beta+n-x)-lgamma(alpha)-lgamma(beta)-lgamma(alpha+beta+n)

def beta_parameters(n,mu):
#     s=((mu*(1-mu))/float(n))**0.5
 #    a=((1-mu)/float(s)-1/float(mu))*(mu**2)

  #   b=a*(1/float(mu) - 1)
     a=1
     b=1
     return a,b








###############################   CLUSTER SIMILARITY     ####################################

def cluster_similarity(cluster1, cluster2):
    sim=cluster(cluster1,cluster2)-cluster1.lik-cluster2.lik
    return sim














