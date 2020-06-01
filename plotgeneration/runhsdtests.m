
% Run LSIR3 tests in HSD
clear all; close all;
addpath('../')

% Problem dimensions
m = 100;
n = 10;

% Set max iterations
maxit = 20;

% Seed random number generator
rng(1);

% Construct RHS in double precision; normalize 
b = randn(m,1);
bnorm = norm(b,2);
b = b./bnorm;

% Set randsvd mode
mode = 3;

% %HSD tests
rng(1);
A = gallery('randsvd',[m,n],1e3,mode);
opdata_lsir_hsd_3 = lsir3(A,b,0,1,2,maxit);
opdata_glsir_hsd_3 = gmreslsir3(A,b,0,1,2,maxit,1e-6);
% opdata_glsir_bd_hsd_3 = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
save('opdata_lsir_hsd_3', 'opdata_lsir_hsd_3');
save('opdata_glsir_hsd_3', 'opdata_glsir_hsd_3');
% save('opdata_glsir_bd_hsd_3', 'opdata_glsir_bd_hsd_3');

rng(1);
A = gallery('randsvd',[m,n],1e4,mode);
opdata_lsir_hsd_4 = lsir3(A,b,0,1,2,maxit);
opdata_glsir_hsd_4 = gmreslsir3(A,b,0,1,2,maxit,1e-6);
% opdata_glsir_bd_hsd_4 = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
save('opdata_lsir_hsd_4', 'opdata_lsir_hsd_4');
save('opdata_glsir_hsd_4', 'opdata_glsir_hsd_4');
% save('opdata_glsir_bd_hsd_4', 'opdata_glsir_bd_hsd_4');

rng(1);
A = gallery('randsvd',[m,n],1e5,mode);
opdata_lsir_hsd_5 = lsir3(A,b,0,1,2,maxit);
opdata_glsir_hsd_5 = gmreslsir3(A,b,0,1,2,maxit,1e-6);
% opdata_glsir_bd_hsd_5 = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
save('opdata_lsir_hsd_5', 'opdata_lsir_hsd_5');
save('opdata_glsir_hsd_5', 'opdata_glsir_hsd_5');
% save('opdata_glsir_bd_hsd_5', 'opdata_glsir_bd_hsd_5');

rng(1);
A = gallery('randsvd',[m,n],1e6,mode);
opdata_lsir_hsd_6 = lsir3(A,b,0,1,2,maxit);
opdata_glsir_hsd_6 = gmreslsir3(A,b,0,1,2,maxit,1e-6);
% opdata_glsir_bd_hsd_6 = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
save('opdata_lsir_hsd_6', 'opdata_lsir_hsd_6');
save('opdata_glsir_hsd_6', 'opdata_glsir_hsd_6');
% save('opdata_glsir_bd_hsd_6', 'opdata_glsir_bd_hsd_6');

% rng(1);
% A = gallery('randsvd',[m,n],1e9,mode);
% opdata_lsir_hsd_9 = lsir3(A,b,0,1,2,maxit);
% opdata_glsir_hsd_9 = gmreslsir3(A,b,0,1,2,maxit,1e-6);
% opdata_glsir_bd_hsd_9 = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
% save('opdata_lsir_hsd_9', 'opdata_lsir_hsd_9');
% save('opdata_glsir_hsd_9', 'opdata_glsir_hsd_9');
% save('opdata_glsir_bd_hsd_9', 'opdata_glsir_bd_hsd_9');

movefile opdata* data