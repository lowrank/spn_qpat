function [ Hn ] = noisy( H, lvl )
 r =  1 + lvl * (2 * rand(size(H)) - 1);
 Hn = H .* r;
end

