
boundary = @(x, y, v)  (1.0 + sin(4 * pi * x));

sigma = struct(...
    'xs', @derenzo,...
    'xa', @(x)shepp_logan(x, 500, 0.02, 0.02) ...
    );

opt = struct('anisotropy', 0.5, 'angle', 64, ...
    'nodes', [0 0; 1 0; 1 1; 0 1]', 'minArea', 1e-4,...
    'sigma', sigma);

%%
tic;
RTE = rte(opt);
t = toc;

fprintf('Caching time %6.2f seconds\n', t);

RTE.setBoundaryCondition(boundary);

tic;
sigmaT   = opt.sigma.xs(RTE.nodes) + opt.sigma.xa(RTE.nodes); ...
sigmaS   = opt.sigma.xs(RTE.nodes);

RTE.sigmaT = sigmaT;
RTE.sigmaS = sigmaS;

t = toc;
fprintf('Coefficient evaluation time %6.2f seconds\n', t);
% The rte routine gives the adjoint solution.
% U_{-} = u(x, -v);
% Reverse the directions by flipping the bottom half to top.
% U_{+} = u(x, v) as the exact solution.
tic;
Un = RTE.ForwardSolve(); 

t =toc;
fprintf('Solution time %6.2f seconds\n', t);