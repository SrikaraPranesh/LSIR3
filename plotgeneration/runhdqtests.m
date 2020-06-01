% Run LSIR3 tests in HDQ

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

% HDQ tests

rng(1);
A = gallery('randsvd',[m,n],1e3,mode);
opdata_lsir_hdq_3 = lsir3(A,b,0,2,4,maxit);
opdata_glsir_hdq_3 = gmreslsir3(A,b,0,2,4,maxit,1e-12);
save('opdata_lsir_hdq_3', 'opdata_lsir_hdq_3');
save('opdata_glsir_hdq_3', 'opdata_glsir_hdq_3');

keyboard
rng(1);
A = gallery('randsvd',[m,n],1e4,mode);
opdata_lsir_hdq_4 = lsir3(A,b,0,2,4,maxit);
opdata_glsir_hdq_4 = gmreslsir3(A,b,0,2,4,maxit,1e-12);
save('opdata_lsir_hdq_4', 'opdata_lsir_hdq_4');
save('opdata_glsir_hdq_4', 'opdata_glsir_hdq_4');


rng(1);
A = gallery('randsvd',[m,n],1e8,mode);
opdata_lsir_hdq_8 = lsir3(A,b,0,2,4,maxit);
opdata_glsir_hdq_8 = gmreslsir3(A,b,0,2,4,maxit,1e-12);
save('opdata_lsir_hdq_8', 'opdata_lsir_hdq_8');
save('opdata_glsir_hdq_8', 'opdata_glsir_hdq_8');

rng(1);
A = gallery('randsvd',[m,n],1e11,mode);
opdata_lsir_hdq_11 = lsir3(A,b,0,2,4,maxit);
opdata_glsir_hdq_11 = gmreslsir3(A,b,0,2,4,maxit,1e-12);
save('opdata_lsir_hdq_11', 'opdata_lsir_hdq_11');
save('opdata_glsir_hdq_11', 'opdata_glsir_hdq_11');

movefile opdata* data

