# Sequential-agglomerative-Clustering
Performs model-based agglomerative clustering as follows:

1. each user in a separate singleton cluster\
2. first 10,000 data points clustered according to a notion of multiplicative posterior similarity\ 
3. each subsequent data point sequentially added\
4. a decision is made whether to add it to existing cluster or create new cluster
