function [ U] = Ufunc( n )

U = zeros(n, n);
for i = 1:n
    U(i,i) = (-1)^i * sqrt(2.0/pi) * gamma(i+0.5)/gamma(i);
end

end