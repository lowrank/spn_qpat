function [ R] = Rfunc( n)

R = zeros(n, n);
for i = 1:n
    for j =1:n
        R(i, j) = (-1)^(i+j-1)* (gamma(i +0.5) * gamma(j-0.5)/pi/gamma(i)/gamma(j)) * (4*j-3)/(i+j-1)/(2*j-2*i-1);
    end
end


end

