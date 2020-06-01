% Generate plot showing the condition numbers of augmented matrix
% and preconditioned augmented matrix with split preconditioning and 
% left preconditioning using block diagonal factors

addpath('../')

% List of 2-norm condition numbers of A to test
kappas = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10];

% Set system dimensions
m = 100;
n = 10;

% For each tested condition number
for i = 1:numel(kappas)
    
    % Seed random number generator
    rng(1);

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
    [~,R] = house_qr_lp(As,0); % half precision via advanpix
    R = (1/mu)*R*diag(1./diag(D));
    R = R(1:n, 1:n);   % RR is trapezoidal factor
    
    % Construct block diagonal preconditioners in single precision
    P1 = single([sqrt(alpha).*eye(m), zeros(m,n); zeros(n,m),(1/sqrt(alpha)).*single(R')]);
    P2 = single([sqrt(alpha).*eye(m), zeros(m,n); zeros(n,m),(1/sqrt(alpha)).*single(R)]);
    P3 = single([eye(m), zeros(m,n); zeros(n,m), (1/alpha)*single(R')*single(R)]);
    
    % Compute "exact" left preconditioned matrix via Advanpix
    extPCleft = mp(P2,64)\(mp(P1,64)\mp(Aug_A,64));
    
    % Compute "exact" split preconditioned matrix via Advanpix
    extPCsplit = (mp(P1,64)\mp(Aug_A,64))/mp(P2,64);
    
    
    % Record infinity norm condition numbers
    condleft(i) = cond(extPCleft,'inf');    % of left preconditioned matrix
    condsplit(i) = cond(extPCsplit,'inf');  % of split preconditioned matrix
    condA(i) =  double(norm(mp(A,64),'inf')*norm(pinv(mp(A,64)),'inf'));    % of A
    condAugA(i) = cond(Aug_A,'inf');    % of augmented matrix
end

% Generate plot
f = figure();
loglog(condA, condAugA, 'k--','MarkerSize',10,'LineWidth',1);
hold on
loglog(condA, condleft, 'rd-','MarkerSize',10,'LineWidth',1);
loglog(condA, condsplit, 'bo-','MarkerSize',10,'LineWidth',1);

xlabel('$\kappa_\infty(A)$','FontSize',10,'Interpreter','latex')
legend({'$\kappa_\infty(\tilde{A})$','$\kappa_\infty(M_2^{-1} M_1^{-1} \tilde{A})$',...
    '$\kappa_\infty( M_1^{-1} \tilde{A}M_2^{-1})$'},...
    'Interpreter','latex', 'Location','NorthWest')
axis([min(condA), max(condA),1e0,1e25])
set(gca,'FontSize',10)

% Save as fig and as pdf
saveas(f, 'fig/pccondnum', 'fig');
saveas(f, 'fig/pccondnum', 'pdf');