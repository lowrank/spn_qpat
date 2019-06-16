function [ L] = Lfunc( n )
L = zeros(n, n);
for i = 1:n
    L(i,i) = 1.0/(4 * i -3);
end

end
