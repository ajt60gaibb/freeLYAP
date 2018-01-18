freeLYAP
====================
A MATLAB implementation of various matrix equation solvers written by Nick Hale, Alex Townsend, and Heather Wilber.

matrixEquationSolver()
====================
Solves the matrix equation `sum_j A{j} X B{j}^T = F` for any number of terms, where `A` and `B` are cell arrays consisting of `mxm` and `nxn` matrices, respectively. The cost of the algorithm depends on the number of terms in the matrix equation. In particular, for more than `3` terms the cost is typically `O((mn)^3)` operations. 

bartelsStewart()
====================
Solves the generalized Sylvester matrix equation given by 
```
AXB' + CXD' = E
```
where `A`, `C` are `mxm` matrices, `B` and `D` are `nxn` matrices, and `E` is an `mxn` matrix. The algorithm is based on a generalization of Bartels-Stewart algorithm `[1]` and costs `O(m^3+n^3)` operations, see `[2]`.

lyap()
===================
Solves the Lyapunov matrix equation given by `A*X + X*A' + Q = 0`, the Sylvester equation `A*X + X*B + C = 0`, and 
the generalized Lyapunov equation `A*X*E' + E*X*A' + Q = 0`.

References
===================
[1] R. H. Bartels & G. W. Stewart, Solution of the matrix equation
AX +XB = C, Comm. ACM, 15, 820–826, 1972.

[2] J. D. Gardiner, A. J. Laub, J. J. Amato, & C. B. Moler, Solution of the
Sylvester matrix equation AXB^T + CXD^T = E, ACM Transactions on
Mathematical Software (TOMS), 18(2), 223-231, 1992.
