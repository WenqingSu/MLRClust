## Randomized Spectral Clustering for Large-Scale Multi-Layer Networks

**MLRClust** performs bias-adjusted spectral clustering for large-scale directed and undirected multi-layer networks using randomization techniques, including random projection and random sampling. Specifically, the random sampling strategy is first used to sparsify the adjacency matrices of all layers. Then, the random projection strategy is used to accelerate the eigen-decomposition of the sum of squared sparsified adjacency matrices. The communities are ﬁnally obtained via k-means clustering of the eigenvectors. The algorithm not only has low time complexity but also saves storage space.


### Installation
Install MLRClust (PDF manual is already bundled; no need to rebuild vignettes)

devtools::install_github("WenqingSu/MLRClust", build_vignettes = FALSE)
