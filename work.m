%WORK Tester for MINRES based iterative refinement for 
%		least squares problems, based on augmented 
%		matrix.

clear all
close all
rng(1);
mp.Digits(34);

m = 10; n = 5;

b = randn(m,1);
bnorm = norm(b,2);
b = b./bnorm;

A = gallery('randsvd',[m,n],1e2,3);
[opdata,x,iter,gmresits] = minreslsir3(A,b,0,2,4,20,1e-12);
[opdata1,x1,iter1,gmresits1] = minreslsir3_bdiag(A,b,0,2,4,20,1e-12);

xact = double(mp(A)\mp(b));
