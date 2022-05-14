# Sequential agglomerative clustering
## Community detection in interaction networks

Performs model-based agglomerative clustering as follows:

1. each user in a separate singleton cluster<br /> 
2. first 10,000 data points clustered according to a notion of multiplicative posterior similarity<br /> 
3. each subsequent data point sequentially added<br /> 
4. a decision is made whether to add it to existing cluster or create new cluster

Perform clustering by running the following:

```bash
cat toy_data.txt | agglomerative_clustering_sequential.py > results.txt
```
