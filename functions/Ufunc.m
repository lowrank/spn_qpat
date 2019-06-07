function [ R] = Ufunc( n )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

R = zeros(n, n);
for i = 1:n
    R(i,i) = (-1)^i * sqrt(2.0/pi) * gamma(i+0.5)/gamma(i);
end

end