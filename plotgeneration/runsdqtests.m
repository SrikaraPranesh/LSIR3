
% Run LSIR3 tests in HSD

addpath('../')
clear all; close all;
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

%SQD tests

rng(1);
A = gallery('randsvd',[m,n],1e4,mode);
opdata_lsir_sdq_4 = lsir3(A,b,1,2,4,maxit);
opdata_glsir_sdq_4 = gmreslsir3(A,b,1,2,4,maxit,1e-12);
save('opdata_lsir_sdq_4', 'opdata_lsir_sdq_4');
save('opdata_glsir_sdq_4', 'opdata_glsir_sdq_4');

rng(1);
A = gallery('randsvd',[m,n],1e6,mode);
opdata_lsir_sdq_6 = lsir3(A,b,1,2,4,maxit);
opdata_glsir_sdq_6 = gmreslsir3(A,b,1,2,4,maxit,1e-12);
save('opdata_lsir_sdq_6', 'opdata_lsir_sdq_6');
save('opdata_glsir_sdq_6', 'opdata_glsir_sdq_6');

rng(1);
A = gallery('randsvd',[m,n],1e7,mode);
opdata_lsir_sdq_7 = lsir3(A,b,1,2,4,maxit);
opdata_glsir_sdq_7 = gmreslsir3(A,b,1,2,4,maxit,1e-12);
save('opdata_lsir_sdq_7', 'opdata_lsir_sdq_7');
save('opdata_glsir_sdq_7', 'opdata_glsir_sdq_7');


rng(1);
A = gallery('randsvd',[m,n],1e13,mode);
opdata_lsir_sdq_13 = lsir3(A,b,1,2,4,maxit);
opdata_glsir_sdq_13 = gmreslsir3(A,b,1,2,4,maxit,1e-12);
save('opdata_lsir_sdq_13', 'opdata_lsir_sdq_13');
save('opdata_glsir_sdq_13', 'opdata_glsir_sdq_13');

movefile opdata* data

