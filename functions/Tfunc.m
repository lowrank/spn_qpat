function [ T] = Tfunc( n )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

T = zeros(n, n);
for i = 1:n
    for j=1:n
        if i == j
            T(i, j) = (2*i-1);
        elseif j == i+1
            T(i,j) = (2*i);
        else
            T(i,j) = (0);
        end
    end
end
end

