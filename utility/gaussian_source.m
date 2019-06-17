function [ ret ] = gaussian_source( x, center, strength, variance)

if nargin < 3
    strength = 1.0;
    variance = 0.1;
end

ret = strength * (1/sqrt(2*pi) / variance) ...
    * exp( -( (x(1,:) - center(1)).^2 + (x(2,:) - center(2)).^2 ) ...
    / (2 * variance^2) );


end

