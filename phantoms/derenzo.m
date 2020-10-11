function [ ret ] = derenzo(x, scale ,contrast)
    if nargin < 2
        scale = 1.0;
        contrast = 8;
    end
    
    N = 500;

    load('derenzo.mat', 'centroids');
    load('derenzo.mat', 'radius');
    
    xmax = max(centroids(:, 1));
    xmin = min(centroids(:, 1));
    ymax = max(centroids(:, 2));
    ymin = min(centroids(:, 2));

    r = max(xmax - xmin, ymax - ymin) * 1.15;

    centroids(:, 1) = (centroids(:, 1) - xmin*1.2) / r;
    centroids(:, 2) = (centroids(:, 2) - ymin*1.2) / r;

    
    [X,Y] =meshgrid(linspace(0,1, N));
    radius = radius / r;

    ret = scale * ones(size(X));

    for i = 1:size(centroids, 1)
        ret = ret + contrast *radius(i)* ((X - centroids(i, 1)).^2 + (Y - centroids(i, 2)).^2 < 2*radius(i)^2);
    end
    
    
%     ret = imgaussian(ret, 4);
    
    
    X = X(:);
    Y = Y(:);
    F = scatteredInterpolant(X,Y,ret(:));
    ret = F(x(1,:), x(2,:));
    
end


