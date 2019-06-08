function [K] = Kfunc( n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    K = zeros(n, 1);
    for i = 1:n
        K(i) = (-1)^(i-1) * fact2(2*i-1)/fact2(2*i-2) / (2*i-1) /(2*i);
    end

end

function ret = fact2(n)
    ret = prod(n:-2:1);
end
