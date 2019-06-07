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
            assert(isfield(opt, 'order'));
            
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
            
        end
        
        function data = genData(obj)
            % The internal data is from the matrix-vector multiplication.
            % The vector is generated from the SPN system. See
            % "./functions" directory. The solution forms a matrix 
            % representing each mode on each point (node).
            
            
            
            
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
    end
    
end

