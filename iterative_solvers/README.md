freeLYAP/iterative_solvers
====================
A MATLAB implementation of various iterative  solvers based on the alternating direction implicit (ADI) method. These methods are highly efficient  for solving the Sylvester matrix equation `AX-XB= F` 
when

(1) A and B are normal matrices,
(2) linear solves involving A or B can be performed cheaply,
(3) the spectra of A and B are contained in sets E and G, respectively, and E and G are either disjoint real intervals or disjoint disks in the complex plane.

See the `Examples_*.m`  files for various examples and descriptions. 
 
 
adi()
====================
Solves the Sylvester matrix equation `AX-XB= F` using the ADI method [3]. 

fadi()
====================
Solves the Sylvester matrix equation `AX-XB= UV' ` in low rank form [2], i.e., finds `Z`, `D`, and `Y` so that 
`ZDY' \approx X'. 

fiadi()
===================
Solves the Sylvester matrix equation `AX-XB= USV' ` in low rank form when `USV' = F` is an approximate SVD [4]. This function is useful when the rank of `F` is moderate to high and the singular values in `diag(S)' decay rapidly. 

getshifts_adi()
===================
Finds the optimal shift parameters for solving the Sylvester matrix equation `AX-XB= F` with `adi()' or `fadi()` when the spectra of A and B are contained in sets E and G, respectively, and E and G are disjoint real intervals [1,3,4]. These shift parameters also work with `fadi()` for finding a low rank solution.

getshifts_smith()
===================
Finds the optimal shift parameter(s) for solving the Sylvester matrix equation `AX-XB= F` with `adi()` or `fadi()` when the spectra of A and B are contained in sets E and G, respectively, and E and G are disjoint disks in the complex plane [4,5]. 

getshifts_fiadi()
===================
Finds the optimal collection of shift parameters for solving the Sylvester matrix equation `AX-XB= USV' ` with `fiadi()` when the spectra of A and B are contained in sets E and G, respectively, and E and G are disjoint real intervals [4]. 

getshifts_fismith()
===================
Finds the optimal collection of shift parameters for solving the Sylvester matrix equation `AX-XB= USV' ` with `fiadi()`  when the spectra of A and B are contained in sets E and G, respectively, and E and G are disjoint disks in the complex plane [4]. 



References
===================

[1] B. Beckermann and A. Townsend, 
 "On the singular values of matrices with displacement structure."
  SIMAX, 38 (2017): pp. 1227-1248. 

[2] Benner, Peter, Ren-Cang Li, and Ninoslav Truhar. 
 "On the ADI method for Sylvester equations." J. of Comp. and App. Math.
 233.4 (2009): 1035-1045.
 
[3] Lu, An, and Eugene L. Wachspress. 
 "Solution of Lyapunov equations by alternating direction implicit iteration." 
  Comp. & Math. with Appl., 21.9 (1991): pp. 43-58.
  
  [4] Townsend, Alex, and Heather Wilber. "On the singular values of matrices 
 with high displacement rank." arXiv preprint arXiv:1712.05864 (2017).
 
[5] Starke, Gerhard. "Near-circularity for the rational Zolotarev 
 problem in the complex plane." J. of approx. theory 70.1 (1992): 115-130.


