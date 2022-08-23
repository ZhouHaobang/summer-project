# summer-project

These are the codes for numerical experienments in my dissertation. We upload codes for three algorithms.

Three algorithms:
1. Alternating Direction Method of Multipliers (ADMM)
2. Singular Value Thresholding (SVT)
3. A modified Cholesky algorithm 

Data:
1. Create_low_rank_matrix.m: A function to generate low rank matrices.
2. Lenna.tif: A real image for tests.

We implement a function to generate low-rank matrices. We implement ADMM and SVT for recovering the generated low-rank matrices
and reconstruct a real image. Test the modified Cholesky algorithm on an indefinite matrix and make it positive definite (semidefinite).


Codes descriptions:

1. Test_ADMM_SDP.m: Run this file to test the ADMM algorithm for solving matrix completion problems on the generated low-rank matrices.

   Test_ADMM_image.m: Run this file to test the ADMM on a real image "Lenna.tif".

   ADMM_SDP.m: A function to implement the ADMM algorithm for solving the SDP reformulation of matrix completion problems.
   
   
2. Test_SVT.m: Run this file to test the SVT on the generated low-rank matrices.
   
   Test_SVT_image.m: Run this file to test SVT on image reconstruction.

   SVT.m: A function to implement the singular value thresholding algorithm.

   SVS.m: A function embedded in the function "SVT" to shrink singular values.
   

3. modchol_ldlt.m: Run this file to test the modified Cholesky algorithm to make an indefinite matrix positive definite (semidefinite).

Notes:

Please correct the repository location for "Lenna.tif" when you run "Test_ADMM_image.m" and "Test_SVT_image.m".
