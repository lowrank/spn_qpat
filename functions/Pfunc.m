function [ P] = Pfunc( n)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
P = eye(n, n);
for i = 1:n
    P(i,i) = 0.5 * gamma(i + n) * gamma(i-0.5)/gamma(i+n-0.5)/gamma(i);
    for j = 1:(i-1)
        P(i,i) = P(i,i) * (1 - 0.5/j);
    end
    
    for j = 1:(n - i)
        P(i,i) = P(i,i) * (1 +0.5 / j);
    end
    
    P(i,i) = P(i,i) * (2*i-1.5);
% for l=1:(i-1)
%     P(i,i) = P(i, i) * (l+i-1) * (l-i+0.5)/(l-i) / (l+i-1.5);
% end
% 
% for l=(i+1):n
%     P(i,i) = P(i,i) * (l+i-1) * (l-i+0.5)/(l-i) / (l+i-1.5);
% end
% 
% P(i,i) = P(i,i) * (2*i-1) * 0.5;

end


end