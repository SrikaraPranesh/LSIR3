% Generate plot showing the condition numbers of augmented matrix
% and preconditioned augmented matrix with left-preconditioner using 
% QR factors

addpath('../')

% List of 2-norm condition numbers of A to test
kappas = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10];

% Set system dimensions
m = 100;
n = 10;

% Get values of u and uf
fp.format = 'h';
[uf,~,~,xmax,~] = float_params('h');chop([],fp);
u = eps('single');

% For each tested condition number
for i = 1:numel(kappas)
    
    % Seed random number generator
    rng(1);
    
    % Generate randsvd matrix A with desired 2-norm condition number
    % Use mode 3: geometrically distributed singular valu15
    A = gallery('randsvd',[m,n],kappas(i),3);
    
    % Compute optimal scalar alpha
    s = svd(A);
    alpha = 2^(-1/2)*min(s);
    
    % Construct (scaled) augmented matrix
    Aug_A = single([alpha.*eye(m), A; A', zeros(n)]);
    
    
    % Compute QR factorization in half precision
    fp.format = 'h'; chop([],fp);
    [~,~,~,xmax,~] = float_params(fp.format);
    D = diag(1./vecnorm(A));
    mu = 0.1*xmax;
    As = chop(mu*A*D);
    [Q,R] = house_qr_lp(As,0); % half precision via advanpix
    R = (1/mu)*R*diag(1./diag(D));
    R = R(1:n, 1:n);   % RR is trapezoidal factor
    Q1 = Q(:,1:n);

    % Construct left preconditioner in single precision
    P1 = single([alpha.*eye(m), single(Q1*R); single(R'*Q1'), zeros(n)]);
    
    % Compute "exact" preconditioned matrix via Advanpix
    extPC = mp(P1,64)\mp(Aug_A,64);
    
    % Record infinity norm condition numbers
    condpc(i) = cond(extPC,'inf');  % of preconditioned matrix
    condAugA(i) = cond(Aug_A,'inf');    % of augmented matrix
    condA(i) =  double(norm(mp(A,64),'inf')*norm(pinv(mp(A,64)),'inf'));    % of A
    % gmn = (m*n*double(uf))/(1-(m*n*double(uf)));
    % const = 2*m*sqrt(n)*gmn;
    const = double(uf);
    bnd(i) = (1+(const*condA(i)))^2;   % bound on kinf(Aug_A) from analysis
    
end

% Generate plot
f = figure();
loglog(condA, condAugA, 'k--','MarkerSize',10,'LineWidth',1);
hold on
loglog(condA, condpc, 'ro-','MarkerSize',10,'LineWidth',1);
hold on
loglog(condA, bnd, 'b-','MarkerSize',10,'LineWidth',1);
hold on
loglog(condA, (1/u).*ones(numel(kappas),1), 'k:','MarkerSize',10,'LineWidth',1);

xlabel('$\kappa_\infty(A)$','FontSize',10,'Interpreter','latex')
legend({'$\kappa_\infty(\tilde{A})$','$\kappa_\infty( M^{-1} \tilde{A})$',...
    '$(1+u_f \kappa_\infty(A))^2$', '$u^{-1}$'},'Interpreter','latex', 'Location','NorthWest')
axis([min(condA), max(condA),1e0,1e15])
set(gca,'FontSize',10)

% Save as fig and as pdf
saveas(f, 'fig/girpccond', 'fig');
saveas(f, 'fig/girpccond', 'pdf');
