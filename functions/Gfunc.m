function [ P ] = Gfunc( n)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
P = zeros(n, n);
for i = 1:n
    for j =1:n
        P(i, j) =  1/((i-0.25)^2 - (j-0.75)^2);
    end
end


end

