function [V,H] = arnoldi(A,v,m,reo)
%ARNOLDI Compute Arnoldi decomposition

if nargin < 4,
    reo = 0;
end;
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

    for r = 0:reo,
        for k = 1:j,
            h = V(:,k)'*w;
            w = w - V(:,k)*h;
            H(k,j) = H(k,j) + h;
        end
    end;
        
    H(j+1,j) = norm(w);
    V(:,j+1) = w/H(j+1,j);
end;

end