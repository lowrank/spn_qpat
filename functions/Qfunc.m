function [ Q] = Qfunc( n)

Q = zeros(n, n);
for i = 1:n
    Q(i,i) = 0.5 * gamma(i + n) * gamma(i+0.5)/gamma(i+n+0.5)/gamma(i);
    for j = 1:(i-1)
        Q(i,i) = Q(i,i) * (1 + 0.5/j);
    end
    
    for j = 1:(n - i)
        Q(i,i) = Q(i,i) * (1 -0.5 / j);
    end
    
    Q(i,i) = Q(i,i) * (2*i-0.5);
end


end