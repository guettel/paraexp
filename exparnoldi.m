function [appr,m,coeff] = exparnoldi(K,M,v,t,m,sigma)
% Approximate expm(t*inv(M)*K)*v using Arnoldi

if nargin < 5, 
    m = 20;
end;
if nargin < 6,
    sigma = 50;
end;

beta = norm(v);
if beta == 0,
    m = 0; appr = v; return
end;

MK = M - K/sigma;
lusolver2();

%A = @(x) MK\(M*x);
A = @(x) lusolver2(MK,M*x);
[V,H] = arnoldi(A,v,m,0);
coeff = expm(t*sigma*(eye(m) - inv(H(1:m,1:m))))*(beta*eye(m,1));
appr = V(:,1:m)*coeff;
