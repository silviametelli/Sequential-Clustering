# Bayesian (sequential) agglomerative clustering
## Community detection in interaction networks

Performs model-based agglomerative clustering as follows:

1. start with each user in a separate singleton cluster<br /> 
2. each possible cluster merge is evaluated according to a notion of multiplicative posterior similarity<br /> 
3. first 10,000 data points clustered in batch according 
4. each subsequent data point sequentially added<br /> 
5. a decision is made whether to add new data point to existing cluster or create new cluster

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
cat toy_data.txt | agglomerative_clustering_sequential.py > results.txt
```
## Reference

Model-Based Clustering and New Edge Modelling in Large Computer Networks. Metelli, Silvia; Heard, Nick. IEEE Intelligence and Security Informatics Conference, Cybersecurity and Big Data. IEEE.
