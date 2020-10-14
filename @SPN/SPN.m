classdef SPN < handle
    %SPN SOLVER FOR SPN APPROXIMATION OF RADIATIVE TRANSFER
    %   THIS CLASS PROVIDES AN INTERFACE TO SOLVE THE SPN SYSTEM, WHICH IS
    %   AN APPROXIMATION TO THE RADIATIVE TRANSPORT EQUATION. HERE N IS THE
    %   ORDER OF SYSTEM. THE SPN SYSTEM IS SOLVED THROUGH THE TRADITIONAL
    %   FINITE ELEMENT METHOD FRAMEWORK. FOR ACCURARCY PURPOSE, THE MESH
    %   MUST BE SUFFICIENT FINE TO RESOLVE THE HIGHEST ORDER'S "MEAN-FREE
    %   PATH".

    %%% PUBLIC PROPOERTIES
    properties (GetAccess = public, SetAccess = private)
        Order    % ORDER OF THE SYSTEM
        Coeff    % COEFFICIENTS OF THE SYSTEM
        Source   % BOUNDARY SOURCE FUNCTION
        Model    % FINITE ELEMENT SOLVER
        Approx   % APPROXIMATION FLAG
        G        % ANISOTROPY PARAMETER
        
    end
    
    
    properties (GetAccess = public, SetAccess = public)
        Reg      % REGULARIZATION PARAMETER
    end
    
    %%% PRIVATE PROPERTIES (DEBUG MODE ON, TURN ACCESS TO PRIVATE WHEN RELEASE)
    properties (GetAccess = public, SetAccess = private)
        cache
    end
    
    %%% FORWARD FUNCTIONS
    methods
        function obj = SPN(opt)
            % SEVERAL ASSERTIONS ON INPUT
            assert(isfield(opt, 'femm_opt'));
            assert(isfield(opt, 'coeff'));
            assert(isfield(opt, 'source'));
            assert(isfield(opt, 'order') && mod(opt.order,2) == 1);
            assert(isfield(opt, 'approx'));
            assert(isfield(opt, 'g'));
            
            if ~(isfield(opt, 'reg'))
                obj.Reg = 1e-8;
            else
                obj.Reg    = opt.reg;
            end
            
            
            % LOAD APPROXIMATION OPTION, NONZERO MEANS USING APPROXIMATION
            obj.Approx = opt.approx;
            obj.Order  = opt.order;
            obj.G      = opt.g;
            obj.Source = opt.source;
            
            
            % LOAD THE FINITE ELEMENT SOLVER
            obj.Model = femm(opt.femm_opt);
            
            % SHORT NAMES FOR A FEW VARIABLES
            nodes = obj.Model.space.nodes;
            elems = obj.Model.space.elems;
            f_ref = obj.Model.facet.ref;
            
           
            % SET COEFFICIENTS ON NODES.
            obj.Coeff.absorption = opt.coeff.absorption(nodes)';
            obj.Coeff.scattering = opt.coeff.scattering(nodes)';
            obj.Coeff.gruneisen  = opt.coeff.gruneisen(nodes)';
            
            % SET COEFFICIENTS ON QUADRATURE-NODES, GRUNEISEN IS NOT
            % INVOLVED HERE.
            qAbsorption = obj.mapping(obj.Coeff.absorption, elems, f_ref');
            qScattering = obj.mapping(obj.Coeff.scattering, elems, f_ref');
            
            
            
            % THE MATRICES USED UNDER SIMPLIFICATION IN PAPER.
            obj.cache.S  = obj.Model.build('s', 1./qScattering);
            obj.cache.MA = obj.Model.build('m', qAbsorption);
            obj.cache.MS = obj.Model.build('m', qScattering);
            obj.cache.Se = obj.Model.build('s', 1); % for regularization?
            obj.cache.M  = obj.Model.build('m', 1); % for L2 norm
            
            % BOUNDARY MATRIX.
            obj.cache.E = obj.Model.build('e', 1, 'all');
            
            % LOAD PRE-DEFINED MATRICES
            T = Tfunc( (obj.Order+1)/2 );
            R = Rfunc( (obj.Order+1)/2 ); 

            obj.cache.K = inv(T);
            obj.cache.B = R * obj.cache.K;
            obj.cache.k = Kfunc( (obj.Order+1)/2  );
            
            
        end
        
        function [M] = assemblePreCondMatrix(obj)
            % RETURN THE ASSEMBLED MATRIX DIRECTLY. SOLVING THE LINEAR
            % EQUATION WILL BE CHANLLENGING IF THE SYSTEM IS LARGE.
            % PARTICULARLY THE PRECONDITIONER IS NOT AVAILABLE EXPLICITLY.
            
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            g = obj.G;
            
            
            M = sparse(N * L, N * L);
            
            
            if obj.Approx == 1
                % INDEXING MIGHT BE SLOW.
                
                % STEP 1. STIFFNESS (diagonal)
                v = diag( 1./( 4 * (1:L) - 1) ./ (1 - g.^(2*(1:L) - 1)) );
                M = kron(v, obj.cache.S);
                
                % STEP 2. MASS 
                % USE KRON TO DEFINE THE MASS MATRIX 
                
                s_1 = obj.cache.K(1, :);
                M = M + kron(s_1' * s_1, obj.cache.MA);
                
                for n = 2:L
                    s_n = obj.cache.K(n, :);
                    M = M + (1-g^(2 * n - 2)) * (4*n-3) * kron(s_n' * s_n, obj.cache.MS);
                end
                
                % STEP 3. TRACE                
                M = M + kron(obj.cache.B, obj.cache.E);
                
            else
                disp('Not Implemented Error.\n');
            end
            
        end

        function [M] = assemblePreCondMatrix_scattering(obj, x)
            % RETURN THE ASSEMBLED MATRIX DIRECTLY. SOLVING THE LINEAR
            % EQUATION WILL BE CHANLLENGING IF THE SYSTEM IS LARGE.
            % PARTICULARLY THE PRECONDITIONER IS NOT AVAILABLE EXPLICITLY.
            elems = obj.Model.space.elems;
            f_ref = obj.Model.facet.ref;
            q_x = obj.mapping(x, elems, f_ref');
            
            Sx  = obj.Model.build('s', 1./q_x);
            MSx = obj.Model.build('m', q_x);
            
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            g = obj.G;
            
            
            M = sparse(N * L, N * L);
            
            
            if obj.Approx == 1
                % INDEXING MIGHT BE SLOW.
                
                % STEP 1. STIFFNESS (diagonal)
                v = diag( 1./( 4 * (1:L) - 1) ./ (1 - g.^(2*(1:L) - 1)) );
                M = kron(v, Sx);
                
                % STEP 2. MASS 
                % USE KRON TO DEFINE THE MASS MATRIX 
                
                s_1 = obj.cache.K(1, :);
                M = M + kron(s_1' * s_1, obj.cache.MA);
                
                for n = 2:L
                    s_n = obj.cache.K(n, :);
                    M = M + (1-g^(2 * n - 2)) * (4*n-3) * kron(s_n' * s_n, MSx);
                end
                
                % STEP 3. TRACE                
                M = M + kron(obj.cache.B, obj.cache.E);
                
            else
                disp('Not Implemented Error.\n');
            end
            
        end        
        
        function [Y_] = assemble(obj, X_)
            % RETURN THE ASSEMBLED MATRIX APPLIED TO A VECTOR/MATRIX INSTEAD OF 
            % RETURNING THE WHOLE MATRIX.
            
            % THE REASON IS: SAVING THE FULL MATRIX MIGHT BE TOO DIFFICULT
            % FOR A LAPTOP TO WORK IF THERE ARE TOO MANY MODES INVOLVED.
            % THE OTHER REASON IS: FOR THE APPROXIMATED VERSION, THE
            % MATRICES ARE REPEATING THEMSELVES, ONE SHOULD NOT CREATE THE
            % MATRICES EACH TIME.
            
            % RESHAPE THE INPUT?
            
            
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            g = obj.G;
            
            % VALIDATION OF INPUT, DEBUG MODE ONLY. COMMENT OUT IF SPEED IS
            % IMPORTANT.
            
            X = reshape(X_, N, L);
            assert(numel(X) == N * L);
            
            % PREPARE THE OUTPUT
            Y = zeros(N, L);
            
            if obj.Approx == 1
                % USE APPROXIMATION MODEL. THE COEFFICIENTS ARE DECOUPLED.
                
                % STEP 1. STIFFNESS
                for n = 1:L
                    Y(:, n)= Y(:, n) +...
                        obj.cache.S * X(:, n)/ (4 * n - 1) / (1 - g^(2*n-1));
                end
                
                % STEP 2. MASS
                
                % PREPARE, IT IS NOT CHEAP HERE.
                Z1 = obj.cache.MS * X;
                Z2 = obj.cache.MA * X;
                
                for n = 1:L
                    s_n = obj.cache.K(n, :); % ROW VECTOR !!
                    theta = zeros(N, 1);
 
                    % GET THE COMMON VECTORS
                    if n ~= 1
                        for k = 1:L
                            theta = theta + s_n(k) * Z1(:, k);
                        end
                    else
                        for k = 1:L
                            theta = theta + s_n(k) * Z2(:, k);
                        end
                    end
                    
                    % PROPORTIONAL APPLY TO EACH COMPONENT. HERE THE
                    % TRIANGULAR SHAPE OF T IS NOT USED, THE TIMING HERE IS
                    % ABOUT TWICE EXPENSIVE. 
                    
                    for j = 1:L
                        if n~= 1
                            % IT SHOULD BE SIGMA_A + (1 - G^(2N-2)) SIGMA_S
                            Y(:, j) = Y(:,j) +  (4*n-3) * (1 - g^(2*n-2)) * s_n(j) * theta;
                        else
                            Y(:, j) = Y(:,j) +  (4*n-3) * 1 * s_n(j) * theta;
                        end
                    end
                end
                        
                % STEP 3. TRACE
                Z3 = obj.cache.E * X;
                
                for j = 1:L
                    for k =  1:L
                       Y(:, j) = Y(:, j) + obj.cache.B(j, k) * Z3(:, k);
                    end
                end
            else
                
                disp('Not Implemented Error.\n');
                
            end
            
            % RESHAPE?
            Y_ = reshape(Y, N*L, 1);
        end

        function [l_] = load(obj, f)
            % GET LOAD VECTOR FROM F.
            load_f = obj.mapping1D(f, obj.Model.space.edges, obj.Model.edge.ref');
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            
            l = zeros(N, L);
            y = obj.Model.build('g', load_f ,'all');
            for i = 1:L
                l(:, i) = obj.cache.k(i) * y;
            end
            
            l_ = reshape(l, N*L, 1);
        end
        
        function plotSolution(obj, Y)
            % PLOT ALL MODES
            L = (obj.Order + 1) / 2;
            N = size(obj.Model.space.nodes, 2);
            
            C = floor(sqrt(L));
            R = floor(L / C);
            
            if C*R < L
                R = R + 1;
            end

            idx = 0;
            for r = 1:R                    
                for c = 1:C
                    idx = idx + 1;
                    if idx > L
                        break;
                    end
                    subplot(R,C,idx);
                    trisurf(obj.Model.space.elems(1:3,:)', ...
                        obj.Model.space.nodes(1,:), obj.Model.space.nodes(2,:),...
                        Y((idx-1)*N+1:idx*N), ...
                        'EdgeColor', 'None');
                    
                    view(2);
                    colormap jet;
                    colorbar;
                    shading interp;
                end

            end
            

        end
        
        function plotData(obj, H)
            % fix the size
            figure('Renderer', 'painters', 'Position', [10 10 600 480])
            trisurf(obj.Model.space.elems(1:3,:)', ...
                    obj.Model.space.nodes(1,:), obj.Model.space.nodes(2,:),...
                    H, ...
                    'EdgeColor', 'None');
            view(2);
            colormap jet;
            colorbar;
            shading interp;
        end
        
    end
    
    %%% BACKWARD FUNCTIONS FOR SINGLE COEFFICIENT
    methods 
        function sigmaA = SingleRecSigmaA(obj, H, f)
            
            
            % DUE TO LINEARITY, ONE NEEDS TO BUILD THE BILINEAR FORM
            % SEPARATELY.
            M = obj.assembleSingleRecSigmaA();
            
            % Get the sigmaA * s1 * Phi by dividing Gruneisen
            h = H ./ obj.Coeff.gruneisen;
            
            elems = obj.Model.space.elems;
            f_ref = obj.Model.facet.ref;
            qh = obj.mapping(h, elems, f_ref');    
            
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            
            l = zeros(N, L);
            
            y = obj.Model.build('l', qh);
            s1 = obj.cache.K(1, :);
            for i = 1:L
                l(:, i) = s1(i) * y;
            end
           
            l_ = reshape(l, N*L, 1);
    
            total_load = f - l_;
            
            x = M \ total_load;
            x_unpack = reshape(x, N, L);
            h_      = (x_unpack * s1');
            
            sigmaA  = h ./ h_;

            
        end
        
        
        function ret  = SingleRecSigmaS(obj, H, f, x0, option)
            N  = size(obj.Model.space.nodes, 2);
%             L  = (obj.Order + 1 ) / 2;
            
            if nargin < 4
                x0 = ones(N,1);
            end
            
            opts    = struct( 'factr', 1e-20, 'pgtol', 1e-12, 'm', 400, 'x0', x0, 'maxIts', 2e2, 'maxTotalIts', 1e5);
            opts.printEvery = 1;
            tic;
            [ret, ~, hist] = lbfgsb_c(@(x)(obj.scatteringRecon_objfcn(x, f, H, option)),...
                zeros(N, 1), inf * ones(N, 1), opts);

            t = toc;
            fprintf('\n\ntotal time in L-BFGS is %2.3e seconds, using %2.3e iterations.\n\n', t, hist.iterations);
            
            
        end
        
        
        function [f, grad] = scatteringRecon_objfcn(obj, sigma_s, f, H_ast, option)
        % OBJECTIVE FUNCTION FOR SCATTERING RECONSTRUCTION.
            h_ast = H_ast ./ obj.Coeff.gruneisen ./ obj.Coeff.absorption;

            N  = size(obj.Model.space.nodes, 2);
            L  = (obj.Order + 1 ) / 2;
            g  = obj.G;


            % BUILD THE MATRIX M IF MEMORY IS ENOUGH (VALID FOR ORDER LESS THAN 19)
            M  = obj.assemblePreCondMatrix_scattering(sigma_s);

            load = obj.load(f);

            x = M\ load;

            x_unpack = reshape(x, N, L); % N x L
            s1     = obj.cache.K(1, :); % 1 x L
            h      = (x_unpack * s1'); 
            
            if option == 1
            mismatch = (obj.cache.M + obj.cache.S) * (h- h_ast); % N x 1
            else
            mismatch = (obj.cache.M) * (h- h_ast);
            end
%             mismatch =  (h- h_ast); % N x 1
            
            f = 0.5 * (h - h_ast)'* mismatch; % L2 norm
           
            f = f + 0.5 * obj.Reg * sigma_s' * obj.cache.Se  * sigma_s;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %               GRADIENT COMPUTATION
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            grad = zeros(N , 1); % gradient with respect to sigma_s only.
            
            V = M \ (kron(s1', mismatch)); % NL x 1
            V_unpack = reshape(V, N, L);
            
            % GRADIENT FROM STIFFNESS MATRIX
            for l = 1:L
                c = 1/( 4 * l - 1) ./ (1 - g^(2*l- 1));
                
                grad = grad + c * ...
                    obj.Model.adj(V_unpack(:, l), x_unpack(:, l), ...
                    1./ (sigma_s.^2), zeros(N, 1));
                
            end
            
            for l = 2:L
                c = (1-g^(2 * l - 2)) * (4*l-3);
                sl = obj.cache.K(l, :);
                sl_x = x_unpack * sl';
                sl_V = V_unpack * sl';
                
                grad = grad - c * ...
                    obj.Model.adj(sl_V, sl_x, ...
                    zeros(N,1), ones(N, 1));
                
            end
            
            grad = grad + obj.Reg * obj.cache.Se * sigma_s;
            
            grad(unique(obj.Model.space.edges)) = 0;
            
            
        end
        
    end
    
    methods (Access = private)
        function [M] = assembleSingleRecSigmaA(obj)
            % RETURN THE ASSEMBLED MATRIX DIRECTLY. SOLVING THE LINEAR
            % EQUATION WILL BE CHANLLENGING IF THE SYSTEM IS LARGE.
            % PARTICULARLY THE PRECONDITIONER IS NOT AVAILABLE EXPLICITLY.
            
            N = size(obj.Model.space.nodes, 2);
            L = (obj.Order + 1) / 2;
            g = obj.G;
            
            
            M = sparse(N * L, N * L);
            
            
            if obj.Approx == 1
                % INDEXING MIGHT BE SLOW.
                
                % STEP 1. STIFFNESS (diagonal)
                v = diag( 1./( 4 * (1:L) - 1) ./ (1 - g.^(2*(1:L) - 1)) );
                M = kron(v, obj.cache.S);
                
                % STEP 2. MASS 
                % USE KRON TO DEFINE THE MASS MATRIX 
                
                %%% The data part replaced into the load.
%                 s_1 = obj.cache.K(1, :);
%                 M = M + kron(s_1' * s_1, obj.cache.MA);
                
                for n = 2:L
                    s_n = obj.cache.K(n, :);
                    M = M + (1-g^(2 * n - 2)) * (4*n-3) * kron(s_n' * s_n, obj.cache.MS);
                end
                
                % STEP 3. TRACE                
                M = M + kron(obj.cache.B, obj.cache.E);
                
            else
                disp('Not Implemented Error.\n');
            end
        end
        
        
        
        
    end
    
    
    
    methods(Static)
        %%% MAPPING FROM NODES TO QUADRATURE-NODES
        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        
        function [interpolate] = mapping1D(func, edges, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(edges, 2));
            for i = 1:size(edges, 2)
                interpolate(:, i) = trans_ref * func(edges(:, i));
            end
        end
    end
    
end

