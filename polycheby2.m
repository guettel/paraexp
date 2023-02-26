function [f,m] = polycheby2(A,v,tol,m_max,a,b)
% POLYCHEBY2  Matrix exponential via Chebyshev interpolant.
% Approximate expm(A)*v to absolute accuracy tol by Chebyshev expansion.
% A can also be given as a function handle.
% The spectrum of A should be contained in an interval [a,b], or it should 
% be close to this interval [a,b].
% For convenience, this code uses FFT to approximate the Chebyshev 
% coefficients (although these are known explicitly as Bessel functions.)
% It also works if [a,b] is a non-real interval, but much more iterations 
% will be required (linearly in the length of the interval).
% If no interval [a,b] is provided then a Gershgorin estimate is used.
% Stefan Guettel, 2010.

if nargin < 3, tol = 1e-6; end
if nargin < 4, m_max = 4096; end
if nargin < 5,
    if isnumeric(A),
        r = sum(abs(A),2) - abs(diag(A)); % Gershgorin
        a = real(full(min(diag(A) - r)));
        b = real(full(max(diag(A) + r)));
    else
        error('A is not numeric, must provide expansion interval [a,b].');
    end
end

% Cheby expansion on [a,b]
g = @(x) exp(((b-a)*x + a + b)/2);
N = 2*m_max;
t = cos( pi*(0:N)/N );
gt = g(t);
c = [ gt , gt(N:-1:2) ];
c = fft(c)/N;
if isreal([a,b]),
    c = real(c);
end
c(1) = c(1)/2; 
%c(m_max) = c(m_max)/2;
c = c(1:m_max);
err = cumsum(abs(c(end:-1:1)));
ind = find(err*norm(v) >= tol,1,'first');
if ind <= 20, % use tail of 20 coefficients for error est.
    error('Maximal number of predicted Chebyshev iterations exceeded, increase m_max in polycheby');
end
if isempty(ind),
    m = 1;
else
    m = m_max - ind;
end

% transform A to [-1,1]
if isnumeric(A),
    AA = @(v) ((2/(b-a))*A*v - ((a+b)/(b-a))*v);
else
    AA = @(v) ((2/(b-a))*A(v) - ((a+b)/(b-a))*v);
end

% Chebyshev recursion
T0 = v;
T1 = AA(v);
f = c(1)*T0;
f = f + c(2)*T1;
for k = 3:m,
    [ T1 , T0 ] = deal( 2*AA(T1) - T0 , T1 );
    f = f + c(k)*T1;
end


