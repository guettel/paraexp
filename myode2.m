function ode = myode2(k,f)
% Defines linear ODE u'' = A*u + g(t) for
% wave equation with oscillating hat source.
% The propagation speed is k and the
% frequency is f. For higher frequencies
% one should also increase N.

N = 100;
ode.A = [ sparse(N,N) speye(N) ; -1*k*(N+1)^2*gallery('tridiag',N) sparse(N,N) ];
x = (1:N)'/(N+1);
w = 0.05; h = 100*sqrt(k);
ode.g = @(t,u) [ zeros(N,1) ; 1*h*max(0, 1-abs((.5 + (.5-w)*sin(2*pi*f*t)) - x)/w)];
ode.u0 = 0*[ .0*ones(N,1) ; 4*x.*(1-x) ];
ode.t = [0,1];

if 0,
    % Neumann
    ode.A(1,1) = ode.A(1,1)/2; 
    ode.A(N,N) = ode.A(N,N)/2;
    ode.u0 = 0*x + 1;
end;
    


