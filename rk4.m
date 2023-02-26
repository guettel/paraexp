function [ts,yout,n] = rk4(func,ts,u0,h)
%[ts,yout] = RK4(func,ts,u0,h)
%   Runge-Kutta 4 time stepping for the function
%   func(t,u), time points ts, initial value u0
%   and step size h.


% initialization
n = 0;
ts = ts(:);
u0 = u0(:);
yout = zeros(length(ts),length(u0));
yout(1,:) = u0;

% coarse time stepping
for j = 1:length(ts)-1,
    stps = ceil((ts(j+1)-ts(j))/h);
    hh = (ts(j+1) - ts(j)) / stps;
    tt = ts(j);
    
    % fine time stepping
    for k = 1:stps,
        k1 = func(tt,u0);
        k2 = func(tt+.5*hh,u0+.5*hh*k1);
        k3 = func(tt+.5*hh,u0+.5*hh*k2);
        k4 = func(tt+hh,u0+hh*k3);
        u0 = u0 + (1/6*hh)*(k1+k4+2*(k2+k3));
        tt = tt + hh;
        n = n + 1;
    end;
    yout(j+1,:) = u0;
    
end;

