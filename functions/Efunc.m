function [ P, Q] = Efunc( n )
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here

P =  eye(n, n);
Q = P;
for i = 1:n
    P(i,i) = sqrt(i + n -0.75) * sqrt(i + n + 0.25) * sqrt(i-0.75) * sqrt(i+0.25) * (i-0.25)/(i-0.75) / i / (i+n);
    Q(i,i) = sqrt(i-0.25) * sqrt(0.75 + n- i)/(0.5 + n -i);
end 

end

