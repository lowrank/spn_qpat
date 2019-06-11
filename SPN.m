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
    
    %%% PRIVATE PROPERTIES (DEBUG MODE ON, TURN ACCESS TO PRIVATE WHEN RELEASE)
    properties (GetAccess = public, SetAccess = private)
        cache
    end
    
    methods
        function obj = SPN(opt)
            % SEVERAL ASSERTIONS ON INPUT
            assert(isfield(opt, 'femm_opt'));
            assert(isfield(opt, 'coeff'));
            assert(isfield(opt, 'source'));
            assert(isfield(opt, 'order') && mod(opt.order,2) == 1);
            assert(isfield(opt, 'approx'));
            assert(isfield(opt, 'g'));
            
            
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
                for n = 1:L
                    M((n-1)*N+1:n*N, (n-1)*N+1:n*N) = ...
                        obj.cache.S / (4 * n - 1) / (1 - g^(2*n-1));
                end
                
                % STEP 2. MASS 
                % SKIPPED SINCE IT IS USUALLY SMALL.
                
                % STEP 3. TRACE                
                for j = 1:L
                    for k =  1:L
                       M((j-1)*N+1:j*N, (k-1)*N+1:k*N) = ...
                           M((j-1)*N+1:j*N, (k-1)*N+1:k*N) + obj.cache.B(j, k) * obj.cache.E;
                    end
                end
                
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
        
        
        function plot(obj, Y)
            % PLOT ALL MODES
            L = (obj.Order + 1) / 2;
            N = size(obj.Model.space.nodes, 2);
            
            C = floor(sqrt(L));
            S = L - C^2;
            
            if S > 0
                R = C + 1;
            else
                R = C;
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

