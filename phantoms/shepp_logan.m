function [ ret ] = shepp_logan( x, N, scale, shift)

%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
if nargin < 2
    N = 500;
    scale = 0.1;
    shift = 0.05;
end

[X,Y] =meshgrid(linspace(0,1, N));
X = X(:);
Y = Y(:);
ph = phantom('Modified Shepp-Logan', N);

% ph = imgaussian(ph, 20);

F = scatteredInterpolant(X,Y,ph(:));
ret = F(x(1,:), x(2,:)) * scale + shift;


end

