% DEMO FOR SPN CLASS

absorption = @(x)(0.5 + 0.1 * x(1,:));
scattering = @(x)(1.5  + 2.2 * x(2,:));
gruneisen  = @(x)(0.75 + 0 * x(1,:));
source     = @(x)(0.2 * x(1,:) .* x(2,:));

%%% FEM CLASS OPTION
femm_opt   = struct(...
    'deg', 1,...
    'qdeg', 3, ...
    'min_area', 4e-5, ...
    'edge', [0 1 1 0; 0 0 1 1]...
    );

%%% coefficients set.
coeff_opt  = struct(...
    'absorption', absorption, ...
    'scattering', scattering, ...
    'gruneisen', gruneisen);

%%% The option struct for SPN.
opt = struct(...
    'order', 1, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source,...
    'approx', 1, ...
    'g', 0.5 ...
    );

Sp = SPN(opt);
N  = size(Sp.Model.space.nodes, 2);
L  = (Sp.Order + 1 ) / 2;

% BUILD THE MATRIX M IF MEMORY IS ENOUGH.
M  = Sp.assemblePreCondMatrix();

f = Sp.Source(Sp.Model.space.nodes)';

load = Sp.load(f);

% tic;
% x = gmres(@Sp.assemble, load, 10, 1e-6, 3200, M, []);
% toc;

% NAIVE UMFPACK SOLVER.
x = M \ load;
Sp.plot(x);

%%% CALCULATE DATA FOR A SINGLE INSTANCE.

x_unpack = reshape(x, N, L);
s1     = Sp.cache.K(1, :);

H      = Sp.Coeff.absorption.* (x_unpack * s1'); % GRUNEISEN ALSO MULTIPLIED.
Sp.plotData(H);

%%% BACKWARD SOLVER.
%
% ADJOINT STATE REQUIRES THE TRANSPOSE OF M, BUT M IS SYMMETRIC (OPERATOR
% IS ADJOINT).



