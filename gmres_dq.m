function [x, error, its, flag] = gmres_dq( A, x, b, L1, L2, R1, R2, restrt, max_it, tol)

% Solves Ax = b by solving the preconditioned linear system 
% (L1*L2)^{-1}A(R1*R2)^{-1}x=(L1*L2)^{-1}b
% using the Generalized Minimal residual ( GMRES ) method.
% Currently uses (preconditioned) residual norm to check for convergence 
% (same as Matlab GMRES)
% Double precision used throughout, except in applying preconditioned matrix 
% to a vector, which is done in quad precision via Advanpix
%
% input   A        REAL nonsymmetric positive definite matrix
%         x        REAL initial guess vector
%         b        REAL right hand side vector
%         L1       REAL first left preconditioner
%         L2       REAL second left preconditioner
%         R1       REAL first right preconditioner
%         R2       REAL second right preconditioner
%         restrt   INTEGER number of iterations between restarts
%         max_it   INTEGER maximum number of iterations
%         tol      REAL error tolerance
%
% output  x        REAL solution vector
%         error    REAL error norm
%         its     INTEGER number of (inner) iterations performed
%         flag     INTEGER: 0 = solution found to tolerance
%                           1 = no convergence given max_it

flag = 0;
its = 0;

% Ensure double working precision
A = double(A);
b = double(b);
x = double(x);

% Cast preconditioners to working precision
L1 = double(L1);
L2 = double(L2);
R1 = double(R1);
R2 = double(R2);

% Compute initial residual 
rtmp = b-A*x;
r = mp(L1,34)\mp(rtmp,34);
r = mp(L2,34)\r;
r = double(r);

bnrm2 = norm(r);
if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end


error(1) = norm( r ) / bnrm2;
if ( error(1) < tol ) return, end

[n,~] = size(A);                                  % initialize workspace
m = restrt;
V(1:n,1:m+1) = zeros(n,m+1);
H(1:m+1,1:m) = zeros(m+1,m);
cs(1:m) = zeros(m,1);
sn(1:m) = zeros(m,1);
e1    = zeros(n,1);
e1(1) = 1.0;

for iter = 1:max_it,                              % begin iteration
    rtmp = b-A*x;
    
    % Apply left preconditioners to vector in quad precision
    r = mp(L2,34)\(mp(L1,34)\mp(rtmp,34));
    % Store result in double precision
    r = double(r);
    
    V(:,1) = r / norm( r );
    s = norm( r )*e1;
    for i = 1:m,                     % construct orthonormal basis via GS
        its = its+1;
        vcur = V(:,i);
        
        % Apply right preconditioners to vector in quad precision
        vcur = mp(R2,34)\(mp(R1,34)\mp(vcur,34));
        
        % Apply matrix and left preconditioners to vector in quad precision
        vcur = mp(L2,34)\(mp(L1,34)\(mp(A,34)*mp(vcur,34)));

        % Store result in double precision
        w = double(vcur);

        for k = 1:i,
            H(k,i)= w'*V(:,k);
            w = w - H(k,i)*V(:,k);
        end
        H(i+1,i) = norm( w );
        V(:,i+1) = w / H(i+1,i);
        for k = 1:i-1,                              % apply Givens rotation
            temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i);
            H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i);
            H(k,i)   = temp;
        end
        [cs(i),sn(i)] = rotmat( H(i,i), H(i+1,i) ); % form i-th rotation matrix
        temp   = cs(i)*s(i);                        % approximate residual norm
        s(i+1) = -sn(i)*s(i);
        s(i)   = temp;
        H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i);
        H(i+1,i) = 0.0;
        error((iter-1)*m+i+1)  = abs(s(i+1)) / bnrm2;
        if ( error((iter-1)*m+i+1) <= tol ),         % update approximation
            y = H(1:i,1:i) \ s(1:i);                 
            addvec = V(:,1:i)*y;
            
            % Apply right preconditioners to vector in quad precision
            addvec = mp(R2,34)\(mp(R1,34)\mp(addvec,34));
            % Store result in double precision
            addvec = double(addvec);                 % and exit
            
            x = x + addvec;
            break;
        end
    end
    
    if ( error(end) <= tol ), break, end
    y = H(1:m,1:m) \ s(1:m);
    addvec = V(:,1:m)*y;
    
    % Apply right preconditioners to vector in quad precision
    addvec = mp(R2,34)\(mp(R1,34)\mp(addvec,34));
    % Store result in double precision
    addvec = double(addvec);
    
    x = x + addvec;                            % update approximation
    rtmp = b-A*x;
  
    % Apply left preconditioners to vector in quad precision
    r = mp(L2,34)\(mp(L1,34)\mp(rtmp,34));           % compute residual
    % Store result in double precision
    r = double(r);

    s(i+1) = norm(r);
    error = [error, s(i+1) / bnrm2];    
    if ( error(end) <= tol ), break, end;           % check convergence
end

if ( error(end) > tol ) flag = 1; end;                 % converged





function [ c, s ] = rotmat( a, b )

%
% Compute the Givens rotation matrix parameters for a and b.
%
if ( b == 0.0 ),
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) ),
    temp = a / b;
    s = 1.0 / sqrt( 1.0 + temp^2 );
    c = temp * s;
else
    temp = b / a;
    c = 1.0 / sqrt( 1.0 + temp^2 );
    s = temp * c;
end