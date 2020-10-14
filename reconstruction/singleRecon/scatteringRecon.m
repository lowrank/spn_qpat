% SINGLE RECONSTRUCTION OF SCATTERING

%% LOAD COEFFICIENT SET 1
absorption = @(x)shepp_logan(x, 500, 0.02, 0.02);
scattering = @derenzo;
gruneisen  = @(x)(0.5 + 0.25 * sin(2 * pi * x(1,:)) .* sin(2 * pi * x(2,:)));

%% SOURCE FUNCTIONS: CONSTANT, POINT GAUSSIAN SOURCES
source1     = @(x)(gaussian_source(x, [0,0], 1.0, 0.1) + ...
                   gaussian_source(x, [1,0], 1.0, 0.1) + ...
                   gaussian_source(x, [1,1], 1.0, 0.1) + ...
                   gaussian_source(x, [0,1], 1.0, 0.1 )); % SOURCE ONLY RESIDES ON THE BOUNDARY
               
source2     = @(x)(gaussian_source(x, [0,0])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source3     = @(x)(gaussian_source(x, [1,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source4     = @(x)(gaussian_source(x, [0,1])); % SOURCE ONLY RESIDES ON THE BOUNDARY
source5     = @(x)(x(1,:) + 1);
% source5     = @(x)(1 .* sin(4 * pi * x(1,:)) + 1);

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

%%% The option struct for SPN.

MAX_ORDER = 17;
noise_lvl = 0.05/sqrt(N);


for f_order = 1:2:MAX_ORDER
opt = struct(...
    'order', f_order, ...
    'femm_opt', femm_opt, ...
    'coeff', coeff_opt,...
    'source', source1 ,... % not used.
    'approx', 1, ...
    'g', 0.5 ...
    );

%% BUILD FORWARD MODEL
ForwardModel = SPN(opt);
N  = size(ForwardModel.Model.space.nodes, 2);
L  = (ForwardModel.Order + 1 ) / 2;


% BUILD THE MATRIX M IF MEMORY IS ENOUGH (VALID FOR ORDER LESS THAN 19)
tic;
M  = ForwardModel.assemblePreCondMatrix();
t = toc;

f = source5(ForwardModel.Model.space.nodes)';
load = ForwardModel.load(f);

tic;
x = M\ load;
t = toc;

%% CALCULATE DATA FOR A SINGLE INSTANCE.
x_unpack = reshape(x, N, L);
s1     = ForwardModel.cache.K(1, :);
H      = ForwardModel.Coeff.gruneisen .* ...
    ForwardModel.Coeff.absorption.* (x_unpack * s1'); % GRUNEISEN ALSO MULTIPLIED.

%%% Plot the data H.
%     ForwardModel.plotData(H);

H = noisy(H, noise_lvl);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%                     RECONSTRUCTION BEGINS
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for b_order = 1:2:MAX_ORDER
        back_opt = struct(...
            'order', b_order, ...
            'femm_opt', femm_opt, ...
            'coeff', coeff_opt,...
            'source', source1,... % not used.
            'approx', 1, ...
            'g', 0.5 ...
            );

        BackwardModel = SPN(back_opt);

        
        BackwardModel.Reg = 1e-8 *N ;
        ret = BackwardModel.SingleRecSigmaS(H, f, ones(N,1), 1);
        err = norm(ret - ForwardModel.Coeff.scattering)/norm(ForwardModel.Coeff.scattering);
        
        fid = fopen('data.txt', 'a+');
        fprintf(fid, '%8d %8d %6.2e\n', f_order, b_order, err);
        fclose(fid);
        
    end
end



