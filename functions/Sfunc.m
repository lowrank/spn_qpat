function [ P] = Sfunc( n)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
P = pi/2* eye(n, n);
for i = 1:n
    for j = 1:(i-1)
        P(i,i) = P(i,i) * (1 + 0.5/j);
    end
    
    for j = 1:(n - i)
        P(i,i) = P(i,i) * (1 -0.5 / j);
    end
end 

end