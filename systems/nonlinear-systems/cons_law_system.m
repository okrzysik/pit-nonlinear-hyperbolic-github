% Abstract class implementing all things related to a system of nonlinear
% hyperbolic conservation laws. 
%
% Where spatial vectors are parsed into and out of functions they are
% expected to be blocked by variable.
%
% Spatial boundaries can be either periodic or constant. Set in pde_pa.bcs
% as either "periodic" or "constant"
%
% Numerical flux can be either Roe or local Lax--Friedrichs. Set in
% disc_pa.num_flux_id as either "ROE" or "LLF"
%
classdef cons_law_system < handle
   
    properties
       mesh_pa  % Parameter structs
       disc_pa
       pde_pa
       
       m % Number of variables in system
       
       linearization_data = []; % Data for linearizations
    end
    
   
    methods
        
    %% Constructor
    function obj = cons_law_system(pde_pa_, m_)        
        obj.pde_pa = pde_pa_;
        obj.m      = m_;
    end
    
    function set.disc_pa(obj, disc_pa_)
        obj.disc_pa  = disc_pa_;
        assert(strcmp(obj.disc_pa.num_flux_id, 'ROE') || strcmp(obj.disc_pa.num_flux_id, 'LLF'), ...
            'disc.num_flux_id = %s not recognised.', obj.disc_pa.num_flux_id)
        
        if strcmp(obj.disc_pa.num_flux_id, 'ROE') 
            assert(isfield(obj.disc_pa, 'delta_smoothing_parameter'), ...
                'ROE flux requires delta_smoothing_parameter for Harten entropy fix must be specified for ROE flux')
        end
    end
    
    function set.mesh_pa(obj, mesh_pa_)
        obj.mesh_pa  = mesh_pa_;
    end
    
    function set.linearization_data(obj, linearization_data_)
        assert(isa(linearization_data_, 'linearization_data'), 'Linearization_data required to be instance of linearization_data')
        
        obj.linearization_data = linearization_data_;
    end
    
    %% Wrapper functions 
    % Forward Euler step
    function q1 = step(obj, t0_idx, q0)
        if ~isempty(obj.linearization_data); obj.linearization_data.t0_idx = t0_idx; end
        
        L  = obj.spatial_discretization(t0_idx, q0);
        q1 = q0 + obj.mesh_pa.dt * L;
    end
    
    % linearized forward Euler step
    function e1 = step_linearized(obj, t0_idx, e0)
        
        assert(~isempty(obj.linearization_data), 'Taking a linearized step requires non-empty linearization_data');
        obj.linearization_data.t0_idx = t0_idx;
        
        L  = obj.spatial_discretization_linearized(t0_idx, e0);
        e1 = e0 + obj.mesh_pa.dt * L;
    end

    % linearized forward Euler step for characteristic variable k
    function w1 = step_linearized_characteristic_variable(obj, t0_idx, w0, k)

        assert(~isempty(obj.linearization_data), 'Taking a linearized step requires non-empty linearization_data');
        
        if t0_idx == obj.linearization_data.nt
            t0_idx = t0_idx-1;
        end

        obj.linearization_data.t0_idx = t0_idx;

        L = obj.spatial_discretization_linearized_characteristic_variable(t0_idx, w0, k);
        w1 = w0 + obj.mesh_pa.dt * L;
        if norm(w1, inf) > 10*norm(w0, inf)
            error('Unstable...')
        end
    end

    
    
    % Takes the absolute value of lambda if |lambda| > delta, but if 
    % smaller it replaces it with a smoothed version ensuring it's bigger 
    % than delta/2. See LeVeque eq. (15.53)
    function phi_delta = smoothed_absolute_value(obj, lambda)
        delta = obj.disc_pa.delta_smoothing_parameter;
        phi_delta = abs(lambda);
        I = find( phi_delta < delta ) ;
        phi_delta(I) = (lambda(I).^2 + delta^2) / (2*delta);
    end
    
    
    %% Functions that must be implemented by child class
    function q0 = initial_condition(obj, x)
       error('Child class must implement "initial_condition"') 
    end
    
    function L = spatial_discretization(obj, t0idx, q0)
       error('Child class must implement "spatial_discretization"') 
    end
    
    function L = spatial_discretization_linearized(obj, t0idx, e0)
       error('Child class must implement "spatial_discretization_linearized"') 
    end

    % k is the index of the characteristic variable to discretize
    function L = spatial_discretization_linearized_characteristic_variable(obj, t0_idx, w0, k)
        error('Child class must implement "spatial_discretization_linearized"') 
    end
    
    % Map from characteristic variables b to primitive variables a. This
    % is achieved by applying R to b.
    % q is the state used to evaluate the eigenvectors.
    function a = right_eigenvector_map(obj, q, b)
        error('Child class must implement "right_eigenvector_map"')
    end
    
    % Map from primitive variables a to characteristic variables b. This
    % is achieved by applying R^-1 to a.
    % q is the state used to evaluate the eigenvectors.
    function b = left_eigenvector_map(obj, q, a) 
        error('Child class must implement "left_eigenvector_map"')
    end

    % Return the kth linearized wave-speed based on current linearization 
    % data. The wave-speeed will be based on the current t0_idx property of 
    % the linearization data.
    % "neg" means from below an interface, i.e., the reconstruction in the 
    % RHS of the cell to the interface's left.
    % "pos" means from above an interface, i.e., the reconstruction in the 
    % LHS of the cell to the interface's right.
    function [lambda_k_neg, lambda_k_pos] = linearized_wave_speed(obj, k) 
        error('Child class must implement "left_eigenvector_map"')
    end
    
    end
    % End of methods

end