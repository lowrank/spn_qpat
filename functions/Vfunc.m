function [ V] = Vfunc( n )
V = zeros(n, n);
for i = 1:n
   V(i,i) = (-1)^i * sqrt(2.0/pi) * gamma(i-0.5)/gamma(i) * (i-0.75);
end

end

