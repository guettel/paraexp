function [ f,m,err ] = siexpmv( A,v,tol,solver,exact )
%RLEXPMV Shift and invert Lanczos for expm(A)*v
%   A should be symmetric!

if nargin < 5,
    exact = 0*v;
end;

if nargin < 4,
    solver = @backslash;
end;

if nargin < 3,
    tol = 0;
end;

% table of optimal poles, together with attained accuracies
poleAccuracy = [
  5.9330e-1  5.0013e-1 ; 5.7314e-1  7.9481e-2 ; 2.0836e+0  2.1864e-2 ;
  3.7720e+0  8.1517e-3 ; 5.6461e+0  3.4616e-3 ; 3.6439e+0  1.2601e-3 ;
  5.2689e+0  4.6288e-4 ; 6.9078e+0  1.8510e-4 ; 5.2689e+0  7.6655e-5 ;
  6.7894e+0  2.8645e-5 ; 8.4028e+0  1.1231e-5 ; 1.0104e+1  4.5924e-6 ;
  8.4028e+0  1.8381e-6 ; 9.9311e+0  7.1477e-7 ; 1.1603e+1  2.8832e-7 ;
  1.3247e+1  1.1891e-7 ; 1.1536e+1  4.6902e-8 ; 1.3095e+1  1.8660e-8 ;
  1.4695e+1  7.6672e-9 ; 1.3095e+1  3.1566e-9 ; 1.4695e+1  1.2449e-9 ;
  1.6301e+1  5.0442e-10 ; 1.7875e+1  2.0625e-10 ; 1.6207e+1  8.4668e-11 ;
  1.7773e+1  3.3855e-11 ; 1.9377e+1  1.3764e-11 ; 2.1005e+1  5.7143e-12 ;
  1.9377e+1  2.2695e-12 ; 2.1005e+1  9.2395e-13 ; 2.2509e+1  3.8462e-13 ;
  2.1005e+1  1.5596e-13 ; 2.2535e+1  6.2356e-14 ; 2.4120e+1  2.5389e-14 ;
  2.5728e+1  1.0492e-14 ; 2.4120e+1  4.2902e-15 ; 2.5699e+1  1.7434e-15 ;
  2.7286e+1  0 ];

% find pole 
m = find(poleAccuracy(:,2)<=tol,1);
xi = poleAccuracy(m,1);
m = 9; xi = 5.3;

%g = @(x) exp(xi*(x-1)./(x+1+eps));
gm = @(X) expm(xi*((X+eye(size(X)))\(X-eye(size(X)))));


% transform A to [-1,1]
M = xi*speye(size(v,1)) - A;
AA = @(v) solver(M , xi*v + A*v);
f = lanczos(gm,AA,v,m,tol);
err = norm(exact-f);

end




function x = backslash(S,b)
% You may use another linear system solver for S*x = b.
% Note that S is symmetric and remains unchanged for all calls of
% this function, so you could/should reuse LU/Cholesky factors.

x = S\b;
end




function g = lanczos(fm,A,v,m,tol)
%LANCZOS Summary of this function goes here
%   Detailed explanation goes here

if isnumeric(A),
    AA = @(v) A*v;
else
    AA = A;
end;

V = zeros(length(v),m+1);
H = zeros(m+1,m);
beta = norm(v);
V(:,1) = v/beta;

for j = 1:m,
    w = AA(V(:,j));
    H(j,j) = V(:,j)'*w;
    w = w - V(:,j)*H(j,j);
    
    if j > 1,
        w = w - V(:,j-1)*H(j,j-1);
        H(j-1,j) = H(j,j-1);
    end;
    
    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j);
end;

f = fm(H(1:j,1:j))*(beta*eye(j,1));
g = V(:,1:j)*f;

end

