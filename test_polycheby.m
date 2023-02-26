

% test polycheby for 1D Laplacian

n = 500;  % nr of interior grid points
t = 0.1;  % time parameter
A = -(n+1)^2*gallery('tridiag',n);
a = -4*(n+1)^2;  % inclusion of spectral interval of A
b = 0;
v = randn(n,1); 
tic
exact = expm(t*A)*v;
toc

% approximate exp(t*A)*v using polycheby
tol = 1e-3;
tic
[f,m] = polycheby(t*A,v,a,b,tol);
toc
norm(f - exact)
