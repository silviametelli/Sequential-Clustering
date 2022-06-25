# Bayesian (sequential) agglomerative clustering
## Community detection in interaction networks

agglomerative_clustering.py performs model-based agglomerative clustering as follows:

1. start with each user in a separate singleton cluster<br /> 
2. each possible cluster merge is evaluated according to a notion of multiplicative posterior similarity<br />
3a decision is made whether to add the data point to an existing cluster or create new cluster


agglomerative_clustering_sequential.py performs model-based agglomerative clustering as follows:

1. perform agglomerative_clustering.py  on the first 10,000 data points (batch)
2. each subsequent data point is sequentially read in from standard input <br /> 
3. a decision is made whether to add the new data point to an existing cluster or create new cluster


## Data

In this example we consider a bipartite graph of users connecting to computers in an enterprise computer network.
Data needs to be in the following format:

```bash
U1	C1,C2,C3,C4
U2	C1,C2,C3,C5
U3	C3,C6
```
## Run clustering algorithm 

Perform clustering and save results by running the following:

```bash
cat Toy_data.txt | agglomerative_clustering_sequential.py > results.txt
```

## Outputs

The algorithm agglomerative_clustering.py returns a list of tuples (User, Cluster membership): the first element of each tuple is the user, and the second the cluster to which the user belongs to.

```bash
[('U1', 1), ('U2', 2), ('U3', 3), ('U4', 3), ('U5', 4), ('U6', 5)]
```

The algorithm agglomerative_clustering_sequential.py returns a txt file ("agglomerative_sequential_output.txt") with associated cluster memberships.

Sample of console printing during the sequential phase:

```bash
Initialising...
[1, 1, 1, 1, 1, 1, 1]
Start sequential...
SEQUENTIAL CLUSTERING: considering line  8
SEQUENTIAL CLUSTERING: considering line  9
SEQUENTIAL CLUSTERING: considering line  10
SEQUENTIAL CLUSTERING: considering line  11
```


## Reference

Model-Based Clustering and New Edge Modelling in Large Computer Networks. Metelli, Silvia; Heard, Nick. IEEE Intelligence and Security Informatics Conference, Cybersecurity and Big Data. IEEE.
