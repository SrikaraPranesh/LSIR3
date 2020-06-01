%WORK1 tests if alpha is actually required.
%      requires chop, float_params, and house_qr_lp to work.

clear all;
close all;

addpath('../')
% Seed random number generator
rng(1);

% List of 2-norm condition numbers of A to test
kappas = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15];

% Set system dimensions
m = 100;
n = 10;

% For each tested condition number
for i = 1:numel(kappas)
    
    % Generate randsvd matrix A with desired 2-norm condition number
    % Use mode 3: geometrically distributed singular values
    A = gallery('randsvd',[m,n],kappas(i),3);
    
    % Compute optimal scalar alpha
    s = svd(A);
    alpha = 2^(-1/2)*min(s);
    
    % Construct (scaled) augmented matrix
    Aug_A = single([alpha.*eye(m), A; A', zeros(n)]);
    
    % Compute QR factorization in half precision
    fp.format = 'h'; chop([],fp);
    [u,~,~,xmax,~] = float_params(fp.format);
    D = diag(1./vecnorm(A));
    mu = 0.1*xmax;
    As = chop(mu*A*D);
    [Q1,R] = house_qr_lp(As,0); % half precision via advanpix
    R = (1/mu)*R*diag(1./diag(D));
    R = R(1:n, 1:n);   % RR is trapezoidal factor
    
    % Construct block diagonal preconditioners in single precision
    M1i = [(1/alpha)*(eye(m)-(Q1*Q1')) Q1*inv(R'); inv(R)*Q1' alpha*inv(R)*inv(R')];
    M2i = [(eye(m)-(Q1*Q1')) Q1*inv(R'); inv(R)*Q1' alpha*inv(R)*inv(R')];
    M3i = [(eye(m)-(Q1*Q1')) Q1*inv(R'); inv(R)*Q1' inv(R)*inv(R')];

    
    PC1 = mp(M1i,64)*mp(Aug_A,64);
    PC2 = mp(M2i,64)*mp(Aug_A,64);
    PC3 = mp(M3i,64)*mp(Aug_A,64);

    
    % Record infinity norm condition numbers
    cnum(i,1) = cond(double(PC1),inf);
    cnum(i,2) = cond(double(PC2),inf);
    cnum(i,3) = cond(double(PC3),inf);
    condA(i) =  double(norm(mp(A,64),'inf')*norm(pinv(mp(A,64)),'inf'));
    condAugA(i) = cond(Aug_A,'inf');  
end

% Generate plot
f = figure();
loglog(condA, condAugA, 'k--','MarkerSize',10,'LineWidth',1);
hold on
loglog(condA, cnum(:,1), 'rx-','MarkerSize',10,'LineWidth',1);
loglog(condA, cnum(:,2), 'bo-','MarkerSize',10,'LineWidth',1);
loglog(condA, cnum(:,3), 'cd-','MarkerSize',10,'LineWidth',1);

xlabel('$\kappa_\infty(A)$','FontSize',16,'Interpreter','latex')
legend({'$\kappa_\infty(\tilde{A})$','Scaled',...
    'Modified scaled','Unscaled'},...
    'Interpreter','latex', 'Location','NorthWest')
axis([min(condA), max(condA),1e0,1e25])
set(gca,'FontSize',16)
