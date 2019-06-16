function [ G ] = Gfunc( n)

G = zeros(n, n);
for i = 1:n
    for j =1:n
        G(i, j) =  1/((i-0.25)^2 - (j-0.75)^2);
    end
end


end

