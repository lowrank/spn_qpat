%% LOAD COEFFICIENT SET 1
absorption = @(x)shepp_logan(x, 500, 0.02, 0.02);
scattering = @derenzo;
gruneisen  = @(x)(0.5 + 0.25 * sin(2 * pi * x(1,:)) .* sin(2 * pi * x(2,:)));

%% SOURCE FUNCTIONS: CONSTANT, POINT GAUSSIAN SOURCES
source1     =  @(x)(gaussian_source(x, [0,0.5], 1.0, 0.2) + ...
                   gaussian_source(x, [0.5,0], 1.0, 0.2) + ...
                   gaussian_source(x, [0.5,1], 1.0, 0.2) + ...
                   gaussian_source(x, [1,0.5], 1.0, 0.2 )); % SOURCE ONLY RESIDES ON THE BOUNDARY
               
source2     = @(x)(gaussian_source(x, [0,0])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source3     = @(x)(gaussian_source(x, [1,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source4     = @(x)(gaussian_source(x, [0,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source5     = @(x)(1 .* sin(4 * pi * x(1,:)) + 1);

%% SETTING OPTIONS
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

phi = [];
runtime = [];

for order = 1:2:35
tic;
%%% The option struct for SPN.
opt = struct(...
    'order', order, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source1 ,...
    'approx', 1, ...
    'g', 0.5 ...
    );

%% BUILD FORWARD MODEL
ForwardModel = SPN(opt);
N  = size(ForwardModel.Model.space.nodes, 2);
L  = (ForwardModel.Order + 1 ) / 2;

% BUILD THE MATRIX M IF MEMORY IS ENOUGH.
% tic;
M  = ForwardModel.assemblePreCondMatrix();
% t = toc;

fprintf('Assemble time %6.2f seconds\n', t);

%%% LOAD VECTOR
f = source5(ForwardModel.Model.space.nodes)';
load = ForwardModel.load(f);
% 
% tic;
% x = gmres(@ForwardModel.assemble, load, 10, 1e-6, 3200, M, []);
% t=toc;

% NAIVE UMFPACK SOLVER.
% tic;
x = M \ load;
% t = toc;

fprintf('Solution time %6.2f seconds\n', t);
% figure(1);
% ForwardModel.plotSolution(x);

x_unpack = reshape(x, N, L);
s1     = ForwardModel.cache.K(1, :);

phi = [phi (x_unpack * s1')]; % size N x 1.

t = toc;

runtime = [runtime t];

end