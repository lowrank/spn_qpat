% This file contains a demo for the class SPN.

absorption = @(x)(0.01 + 0.01 * x(1,:));
scattering = @(x)(1.5  + 2.2 * x(2,:));
gruneisen  = @(x)(0.75 + 0 * x(1,:));
source     = @(x)(1.0  + 0.2 * x(1,:));

%%% finite element method option struct.
femm_opt   = struct(...
    'deg', 1,...
    'qdeg', 3, ...
    'min_area', 1e-4, ...
    'edge', [0 1 1 0; 0 0 1 1]...
    );

%%% coefficients set.
coeff_opt  = struct(...
    'absorption', absorption, ...
    'scattering', scattering, ...
    'gruneisen', gruneisen);

%%% The option struct for SPN.
opt = struct(...
    'order', 11, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source,...
    'approx', 1, ...
    'g', 0.5 ...
    );

Sp = SPN(opt);
N  = size(Sp.Model.space.nodes, 2);
L  = (Sp.Order + 1 ) / 2;
M  = Sp.assemblePreCondMatrix();

f = Sp.Source(Sp.Model.space.nodes)';

load = Sp.load(f);

tic;
x = gmres(@Sp.assemble, load, 10, 1e-6, 3200, M, []);
toc;

Sp.plot(x);

%%% Calculate the forward solution.


%%% Calculate the backward solution.
