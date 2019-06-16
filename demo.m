% DEMO FOR SPN CLASS WITH APPROXIMATION

% load phantom into coefficients

absorption = @shepp_logan;
scattering = @derenzo;
gruneisen  = @(x)(0.75 + 0 * x(1,:));
source     = @(x)(0.2 * x(1,:) .* x(2,:)); % SOURCE ONLY RESIDES ON THE BOUNDARY

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
    'order', 11, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source,...
    'approx', 1, ...
    'g', 0.5 ...
    );

ForwardModel = SPN(opt);
N  = size(ForwardModel.Model.space.nodes, 2);
L  = (ForwardModel.Order + 1 ) / 2;

% BUILD THE MATRIX M IF MEMORY IS ENOUGH.
tic;
M  = ForwardModel.assemblePreCondMatrix();
t = toc;

fprintf('Assemble time %6.2f seconds\n', t);

%%% LOAD VECTOR
f = ForwardModel.Source(ForwardModel.Model.space.nodes)';
load = ForwardModel.load(f);

% tic;
% x = gmres(@Sp.assemble, load, 10, 1e-6, 3200, M, []);
% toc;

% NAIVE UMFPACK SOLVER.
tic;
x = M \ load;
t = toc;

fprintf('Solution time %6.2f seconds\n', t);
ForwardModel.plotSolution(x);

%%% CALCULATE DATA FOR A SINGLE INSTANCE.

x_unpack = reshape(x, N, L);
s1     = ForwardModel.cache.K(1, :);
H      = ForwardModel.Coeff.gruneisen .* ...
    ForwardModel.Coeff.absorption.* (x_unpack * s1'); % GRUNEISEN ALSO MULTIPLIED.

%%% Plot the data H.
ForwardModel.plotData(H);

H = noisy(H, 0.05);

%%% BACKWARD SOLVER.
%
% ADJOINT STATE REQUIRES THE TRANSPOSE OF M, BUT M IS SYMMETRIC (OPERATOR
% IS ADJOINT).

opt = struct(...
    'order', 1, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source,...
    'approx', 1, ...
    'g', 0.5 ...
    );

BackwardModel = SPN(opt);

%%% 1. Single Reconstruction for absorption.
load = BackwardModel.load(f);
sigmaA = BackwardModel.SingleRecSigmaA(H, load);

fprintf('Reconstruction error is %6.2e\n', ...
    (norm(sigmaA - ForwardModel.Coeff.absorption) / norm(ForwardModel.Coeff.absorption)));


