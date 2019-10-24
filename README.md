## Sparse Matrix Vector multiplication with custom Fault Injection in elements of the vector

The code computes y = A * x (A matrix of n x m size and x and y vectors of n size), in every iteration x being replace by the newly computed vector y. The code is based on the Sparse Matrix-Vector Multiplication FPGA benchmark from the OpenDwarfs project (https://github.com/vtsynergy/OpenDwarfs).

### Data structures

The matrix is stored in the CSR format, with one array for the non-zero elements (floating-point) Ax, one array for the cumulated number of elements per row Ap, and one array for the column indices for each element Aj.

The vectors x and y are simple contiguous arrays of floating-point values.

