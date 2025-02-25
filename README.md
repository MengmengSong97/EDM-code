This is for the code on our paper 

"On the local and global minimizers of the smooth stress function in Euclidean Distance Matrix problems"

---------------------------------------------------------------------------------------------------

Code and functions:

main.m: generates Pbar and Phat either random examples or our Example 1 and Example 2, and then calls lngminTR.m to solve f_L(L).

lngminTR.m: TR method for solving 'f_L(L)' to check if it numerically converges to a local nonglobal minimizer (lngm).

lngminFRobjgradHess.m: calculations for optimal value and derivatives of 'f_L(L)'.

GS.m: Schmidt orthogonal process for getting 'V' instead using 'null(e')'

---------------------------------------------------------------------------------------------------

Data:

For Example 1:

Lbars.mat: 'Lbar' for the n=50, d=1 example 

Llngs.mat: 'Llng' for the n=50, d=1 example, the numerical lngm.

For Example 2:

lbars2.mat: 'lbar2' for the n=100, d=2 example 

llngs2.mat: 'llng2' for the n=100, d=2 example, the numerical lngm.

---------------------------------------------------------------------------------------------------

Illustration of Examples with lngm:

Example 1.m: Data calculation and error analysis in Example 1.

Example 2.m: Data calculation and error analysis in Example 2.
