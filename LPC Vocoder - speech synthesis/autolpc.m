function [A,G,a,r]=autolpc(x,p)

% Inputs:
%   x is the signal frame (usually speech weighted by a window)
%   p is the lpc model order

% Outputs:
%   A is the denominator vector for the lpc solution, i.e.,
%   A=1-a1z^{-1}-a2z^{-2}-...-apz^{-p}
%   G is the lpc model gain (rms prediction error)
%   a is the lpc polynomial (without the 1 term)
%   r is the vector of autocorrelation coefficients

x = x(:);
L = length(x);
r = zeros(p+1,1);
for i=0:p
   r(i+1) = x(1:L-i)' * x(1+i:L);
end
R = toeplitz(r(1:p));
a = R\r(2:p+1);    %<--- compute inv(R)
A = [1; -a];
G = sqrt(sum(A.*r));
end