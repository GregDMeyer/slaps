#SLAPS
####Sparse Linear Algebra on a Partitioned global address Space

SLAPS uses [UPC++](https://bitbucket.org/berkeleylab/upcxx/wiki/Home) to implement sparse matrix-dense vector multiplication (SpMV) on distributed memory. Its features include:

 - **Implicit** remote memory reading and writing through the `Vec` class
 - **Efficient** SpMV through one-sided PGAS communication, in some cases outperforming PETSc's implementation
 - **Header-only** C++ library for ease of use and templated API

Read `./docs/writeup.pdf` to learn more!

For an example, see `./benchmark/slaps`.

&copy; Greg Meyer, 2018.