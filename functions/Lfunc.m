function [ T] = Lfunc( n )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

T = zeros(n, n);
for i = 1:n
    T(i,i) = 1.0/(4 * i -3);
end

end
