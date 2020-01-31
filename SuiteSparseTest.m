%SUITESPARSETEST Tests performance of gmreslsq-ir on matrices from
%   SuiteSparse Matrix collection



clear all; close all;
rng(1);
% Input parameters
maxit = 15;
scale.theta = 0.1; scale.pert = 15; % diagonal perturbation constant
fp.format = 'h'; % low precision format to be considered

% & index.nrows > index.ncols
index = ssget;
indlist = find(index.isReal == 1 & ...
    index.nrows >= 20 & index.nrows <= 2000 & ...
    index.ncols < 400 & ...
    index.sprank == index.ncols);
[nlist,i] = sort(index.nrows(indlist)) ;
indlist   = indlist(i);
nn = length(indlist);
eval_ctest = zeros(nn,4);
mn = zeros(nn,1);
nlist = nlist';
ctest = zeros(nn,1);

fid1 = fopen('gmreslsqir_test.txt','w');
[u,xmins,xmin,xmax,p,emins,emin,emax] = float_params(fp.format);
chop([],fp);

a = 0;
for j=1:nn
    fprintf('Processing matrix %d || Total matrices %d\n',j,nn);
    Problem = ssget(indlist(j));
    A = full(Problem.A);
    %
    AbsA = abs(A);
    mel(j,1) = max(max(AbsA));
    mel(j,2) = min(AbsA(AbsA>0));
    
    
    [d,r(1,1)] = size(A);
    r(1,2) = rank(A);

    if ((r(1,1) == r(1,2)) && (d > r(1,1)))
        a = a+1;
        
        act_ind(a,1) = j;
        rows(a,1) = d;
        eval_ctest(a,1) = r(1,1);
        eval_ctest(a,4) = d;
        eval_ctest(a,3) = cond(A);
        
        
        % GMRES-IR Test
        [m,n] = size(A);
        b = randn(m,1);

        % (half,single,double)
        [S_dat1(a),~,its{1,1}(a,1), t_gmres_its{1,1}(a,1)] = gmreslsir3(A,b,0,1,2,maxit,1e-6);
        [S_dat2(a),~,its{1,2}(a,1),t_gmres_its{1,2}(a,1)] = gmreslsir3_bdiag(A,b,0,1,2,maxit,1e-6);
        [~,~,its{1,3}(a,1)] = lsir3(A,b,0,1,2,maxit);

        
        % (half,double,quad)
        [D_dat1(a),~,its{2,1}(a,1),t_gmres_its{2,1}(a,1)] = gmreslsir3(A,b,0,2,4,maxit,1e-12);
        [D_dat2(a),~,its{2,2}(a,1),t_gmres_its{2,2}(a,1)] = gmreslsir3_bdiag(A,b,0,2,4,maxit,1e-12);
        [~,~,its{2,3}(a,1)] = lsir3(A,b,0,2,4,maxit);

    end
end

% print matrix properties
fprintf(fid1,'Properties of test matrices \n');
for j=1:a
    jj = act_ind(j,1);
    mi = indlist(jj);
    fprintf(fid1,'%d & %s &(%d,%d)& %6.2e & %6.2e & %6.2e\\\\\n',...
        j,index.Name{mi,1},eval_ctest(j,4),eval_ctest(j,1),eval_ctest(j,3),...
        mel(j,1),mel(j,2));
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');


% creating a text file to print the GMRES iteration table
fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'half, single, double combination\n');
for i = 1:a
    t1 = its{1,1}(i,1);
    t2 = its{1,2}(i,1);
    t3 = its{1,3}(i,1);
    
    t1a  =  t_gmres_its{1,1}(i,1);
    t2a  =  t_gmres_its{1,2}(i,1);
    fprintf(fid1,'%d & %d & %d &(%d) & %d &(%d)\\\\ \n',i,t3,t1,t1a,t2,t2a);
end
fprintf(fid1,'\n'); fprintf(fid1,'\n');

fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'half, double, quad combination\n');
for i = 1:a
    t1 = its{2,1}(i,1);
    t2 = its{2,2}(i,1);
    t3 = its{2,3}(i,1);
    
    t1a  =  t_gmres_its{2,1}(i,1);
    t2a  =  t_gmres_its{2,2}(i,1);
    fprintf(fid1,'%d & %d & %d &(%d) & %d &(%d)\\\\ \n',i,t3,t1,t1a,t2,t2a);
end


fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'half, single, double Condition numbers\n');
for i = 1:a
    t1 = S_dat1(i).condAugA;
    t2 = S_dat1(i).condMA;
    t3 = S_dat2(i).condMA;

    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e \\\\ \n',i,t1,t2,t3);
end


fprintf(fid1,'\n'); fprintf(fid1,'\n');
fprintf(fid1,'half, double, quad Condition numbers\n');
for i = 1:a
    t1 = D_dat1(i).condAugA;
    t2 = D_dat1(i).condMA;
    t3 = D_dat2(i).condMA;

    fprintf(fid1,'%d & %6.2e & %6.2e & %6.2e \\\\ \n',i,t1,t2,t3);
end

fclose(fid1);
