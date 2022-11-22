# LJKO

Matlab implementation of the variational upwind finite volume scheme for Wasserstein gradient flows presented in:

[1] C. Cancès, T. O. Gallouët, G. Todeschi. [A variational finite volume scheme for Wasserstein gradient flows.](https://arxiv.org/abs/1907.08305)
Numerische Mathematik, 146(3), pp. 437-480, 2020.

The time discretization is based on the renowned JKO (Jordan-Kinderlehrer-Otto) time discretization and on an implicit linearization of the
Wasserstein distance expressed thanks to the Benamou-Brenier formula. The space discretization
relies on upstream mobility two-point flux approximation finite volumes. The variational structure of the
continuous model is preserved at the discrete level by following a first discretize then optimize approach. The use of upstream mobility enables
to solve efficiently the discrete problem thanks to a Newton scheme.

## Definition of the gradient flow
The scheme can compute the gradient flow with respect to any strictly convex energy functional,
which has to be set via a specific function definition.

**For simplicity, the code is designed to consider only gradient flows with respect to one species (i.e., the energy depends only on one measure).
Tackling multiple species is possible (see [1]) but needs to be done ad hoc.**

## Mesh structure
Two types of meshes are available: \
1 -> regular triangulation of the domain, with only acute angles
     (https://www.i2m.univ-amu.fr/fvca5/benchmark/Meshes/index.html) \
2 -> cartesian grids \
For each mesh, five levels of refinement h_i, 1->5, are available. \
Mesh structure: \
nodes -> array of nodes coordinates [x y] \
cells -> array of cells nodes [#nodes node1 node2 node3 ...] \
edges -> array of edges [node1 node2 K L d_sigma m_sigma m_sigma/d_sigma] \
ind -> indices struct: ind.internal indices of internal edges, ind.bound indices of boundary edges \
cc -> array of cell centers coordinates [x y] \
area -> array of cell measures \
mid -> array of midpoints coordinates [x y] \
h -> meshsize, max(diam(K)) or max(sqrt(area)) \
