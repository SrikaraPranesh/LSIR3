%TABLE_RESULTS generates a latex table of results reported in the paper.

clear all;
close all;

% seed random number generator
rng(1);

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

fid1 = fopen('lsqr_randsvd_ir_steps.txt','w');
fid2 = fopen('lsqr_randsvd_condnum.txt','w');

% condition numbers for various precision combinations
k_hdq = [1e2;1e4;1e7;1e9;1e10;1e11;1e12]; % half,double,quad
k_hsd = [1e3;1e4;1e5;1e6;1e7;1e8]; % half,single,double
k_sdq = [1e3;1e5;1e7;1e9;1e11;1e13;1e15;1e16]; % single,double,quad


% half, double, quad combination
fprintf(fid1,'half, double, quad combination \n');
fprintf(fid2,'half, double, quad combination \n');
for i=1:length(k_hdq)
	A = gallery('randsvd',[m,n],k_hdq(i,1),mode);
	[~,~,iter(i,1)] = lsir3(A,b,0,2,4,maxit); 
	% [opdata2,~,iter(i,2),gmresits(i,2)] = ...
	% 	gmreslsir3(A,b,0,2,4,maxit,1e-12,2);
	[opdata3,~,iter(i,3),gmresits(i,3)] = ...
		gmreslsir3_bdiag(A,b,0,2,4,maxit,1e-12);
	[opdata4,~,iter(i,4),gmresits(i,4)] = ...
		gmreslsir3(A,b,0,2,4,maxit,1e-12,1);
	[~,~,iter(i,6),gmresits(i,6)] = ...
		minreslsir3_bdiag(A,b,0,2,4,maxit,1e-12);

	% print the results into a file
	fprintf(fid1,'%6.2e & %d & %d & (%d)& %d & (%d) & %d & (%d)\\\\\n',...
		k_hdq(i,1),iter(i,1),iter(i,4),gmresits(i,4),iter(i,3),gmresits(i,3)...
							,iter(i,6),gmresits(i,6));

	fprintf(fid2,'%6.2e & %6.2e & %6.2e & %6.2e \\\\\n',...
		k_hdq(i,1),opdata4.condAugA,opdata4.condMA,opdata3.condMA);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid2,'\n'); fprintf(fid2,'\n');

% half, single, double combination
fprintf(fid1,'half, single, double combination \n');
fprintf(fid2,'half, single, double combination \n');
for i=1:length(k_hsd)
	A = gallery('randsvd',[m,n],k_hsd(i,1),mode);
	[~,~,iter(i,1)] = lsir3(A,b,0,1,2,maxit); 
	% [opdata2,~,iter(i,2),gmresits(i,2)] = ...
	% 	gmreslsir3(A,b,0,1,2,maxit,1e-6,2);
	[opdata3,~,iter(i,3),gmresits(i,3)] = ...
		gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
	[opdata4,~,iter(i,4),gmresits(i,4)] = ...
		gmreslsir3(A,b,0,1,2,maxit,1e-6,1);
	[~,~,iter(i,6),gmresits(i,6)] = ...
		minreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);

	% print the results into a file
	fprintf(fid1,'%6.2e & %d & %d & (%d)& %d & (%d) & %d & (%d)\\\\\n',...
		k_hsd(i,1),iter(i,1),iter(i,4),gmresits(i,4),iter(i,3),gmresits(i,3)...
							,iter(i,6),gmresits(i,6));

	fprintf(fid2,'%6.2e & %6.2e & %6.2e & %6.2e \\\\\n',...
		k_hsd(i,1),opdata4.condAugA,opdata4.condMA,opdata3.condMA);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid2,'\n'); fprintf(fid2,'\n');

% single, double, quad combination
fprintf(fid1,'single, double, quad combination \n');
fprintf(fid2,'single, double, quad combination \n');
for i=1:length(k_sdq)
	A = gallery('randsvd',[m,n],k_sdq(i,1),mode);
	[~,~,iter(i,1)] = lsir3(A,b,1,2,4,maxit); 
	% [opdata2,~,iter(i,2),gmresits(i,2)] = ...
	% 	gmreslsir3(A,b,1,2,4,maxit,1e-12,2);
	[opdata3,~,iter(i,3),gmresits(i,3)] = ...
		gmreslsir3_bdiag(A,b,1,2,4,maxit,1e-12);
	[opdata4,~,iter(i,4),gmresits(i,4)] = ...
		gmreslsir3(A,b,1,2,4,maxit,1e-12,1);
	[~,~,iter(i,6),gmresits(i,6)] = ...
		minreslsir3_bdiag(A,b,1,2,4,maxit,1e-12);

	% print the results into a file
	fprintf(fid1,'%6.2e & %d & %d & (%d)& %d & (%d) & %d & (%d)\\\\\n',...
		k_sdq(i,1),iter(i,1),iter(i,4),gmresits(i,4),iter(i,3),gmresits(i,3)...
							,iter(i,6),gmresits(i,6));

	fprintf(fid2,'%6.2e & %6.2e & %6.2e & %6.2e \\\\\n',...
		k_sdq(i,1),opdata4.condAugA,opdata4.condMA,opdata3.condMA);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid2,'\n'); fprintf(fid2,'\n');




fclose(fid1);
fclose(fid2);







