% DEMO FOR SPN CLASS WITH APPROXIMATION

% load phantom into coefficients
absorption = @(x)shepp_logan(x, 500, 0.02, 0.02);
scattering = @derenzo;
gruneisen  = @(x)(0.5 + 0.25 * sin(2 * pi * x(1,:)) .* sin(2 * pi * x(2,:)));

%%% POINT GAUSSIAN SOURCES
source1     = @(x)(gaussian_source(x, [0,0], 1.0, 0.1) + ...
                   gaussian_source(x, [1,0], 1.0, 0.1) + ...
                   gaussian_source(x, [1,1], 1.0, 0.1) + ...
                   gaussian_source(x, [0,1], 1.0, 0.1 )); % SOURCE ONLY RESIDES ON THE BOUNDARY
               
source2     = @(x)(gaussian_source(x, [0,0])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source3     = @(x)(gaussian_source(x, [1,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source4     = @(x)(gaussian_source(x, [0,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source5     = @(x)(x(1,:)+ 1);

%%% FEM CLASS OPTION
femm_opt   = struct(...
    'deg', 1,...
    'qdeg', 4, ...
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
    'order',1, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source1 ,... % not used.
    'approx', 1, ...
    'g', 0.5, ...
    'reg', 1e-10...
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
f = source5(ForwardModel.Model.space.nodes)';
load = ForwardModel.load(f);
% 
% tic;
% x = gmres(@ForwardModel.assemble, load, 10, 1e-6, 3200, M, []);
% t=toc;

% NAIVE UMFPACK SOLVER.
tic;
x = M \ load;
t = toc;

fprintf('Solution time %6.2f seconds\n', t);
% figure(1);
% ForwardModel.plotSolution(x);

%%% CALCULATE DATA FOR A SINGLE INSTANCE.7

x_unpack = reshape(x, N, L);
s1     = ForwardModel.cache.K(1, :);
H      = ForwardModel.Coeff.gruneisen .* ...
    ForwardModel.Coeff.absorption.* (x_unpack * s1'); % GRUNEISEN ALSO MULTIPLIED.



%%% Plot the data H.

% ForwardModel.plotData(H);
% noise_lvl = 0.005;
H = noisy(H, 0.05/sqrt(N));
% res = H2 - H;
% scal = sqrt((res' * ( ForwardModel.cache.M + ForwardModel.cache.S )* res ) / ...
%     (H' * ( ForwardModel.cache.M + ForwardModel.cache.S )* H ));
% 
% H = H + (H2 - H) * noise_lvl / scal;

%%
%%% BACKWARD SOLVER.
%
% ADJOINT STATE REQUIRES THE TRANSPOSE OF M, BUT M IS SYMMETRIC (OPERATOR
% IS ADJOINT).

%%
opt = struct(...
    'order',1, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source1,...
    'approx', 1, ...
    'g', 0.5 ...
    );

BackwardModel = SPN(opt);

% %%% 1. Single Reconstruction for absorption.
% load = BackwardModel.load(f);
% sigmaA = BackwardModel.SingleRecSigmaA(H, load);
% 
% fprintf('Reconstruction error is %6.2e\n', ...
%     (norm(sigmaA - ForwardModel.Coeff.absorption) / norm(ForwardModel.Coeff.absorption)));
% 
% % figure(2);
% % ForwardModel.plotData(ForwardModel.Coeff.absorption);
% % 
% % figure(3);
% % BackwardModel.plotData(sigmaA);
% % 
% % figure(4);
% BackwardModel.plotData(sigmaA - ForwardModel.Coeff.absorption);

%%
BackwardModel.Reg =1e-8 *N ;
ret = BackwardModel.SingleRecSigmaS(H, f, ones(N,1), 1);
BackwardModel.plotData(ret)
norm(ret - ForwardModel.Coeff.scattering)/norm(ForwardModel.Coeff.scattering)

