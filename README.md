Codes and Data in the paper "On the local and global minimizers of the smooth stress function in Euclidean Distance Matrix problems", which is posted at https://arxiv.org/abs/2408.07256

CODE:


main.m: generates Pbar and Phat either random examples or our Example 1 and Example 2, and then calls lngminTR.m to solve f_L(L).


lngminTR.m: TR method for solving 'f_L(L)' to check if it numerically converges to a local nonglobal minimizer (lngm).


lngminFRobjgradHess.m: calculations for optimal value and derivatives of 'f_L(L)'.


GS.m: Schmidt orthogonal process for getting 'V' instead of using 'null(e')'



DATA:

For Example 1:

Lbars.mat: 'Lbar' for the n=50, d=1 example 

Llngs.mat: 'Llng' for the n=50, d=1 example, the numerical lngm.

For Example 2:

lbars2.mat: 'lbar2' for the n=100, d=2 example 

llngs2.mat: 'llng2' for the n=100, d=2 example, the numerical lngm.



ILLUSTRATION:

Example 1.m: Data calculation and error analysis in Example 1.

Example 2.m: Data calculation and error analysis in Example 2.
