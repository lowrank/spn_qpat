function [ S] = Sfunc( n)

S = pi/2* eye(n, n);
for i = 1:n
    for j = 1:(i-1)
        S(i,i) = S(i,i) * (1 + 0.5/j);
    end
    
    for j = 1:(n - i)
        S(i,i) = S(i,i) * (1 -0.5 / j);
    end
end 

end