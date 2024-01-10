% Class defining the parameters used to linearize the scalar conservation
% law in cons_law_scalar.m
%
%
classdef scalar_linearization_pa < handle
 
    properties
       
        % Strategy for linearization of WENO weights (if using WENO).
        weno_linearization; 
        % The options are:
        % 'picard'; % Freezes weights 
        % 'newton'; % Uses gradient of weights
        % 'newton-FD-approx'; % Approximates gradients with FD
        fd_weno_epsilon = 1e-1; % FD step size when approximating gradients of WENO weights using a finite difference
        
        % In the linearized problem there is a lower and upper reconstructions of 
        % the wave-speed at a given interface. If true, this option uses an average
        % average of these two in the numerical flux. In the few tests I
        % tried, this really upset the iteration... 
        average_linearized_interface_wave_speed = ~true;
        
    end
    
    % Constructor
    methods 
        
        % Create class only by setting the weno_linearization since this is
        % realistically the only thing the user should be setting.
        function obj = scalar_linearization_pa(weno_linearization_)
            obj.weno_linearization = weno_linearization_;
        end
    
        function set.weno_linearization(obj, weno_linearization_)
            assert(strcmp(weno_linearization_, 'picard') ...
                || strcmp(weno_linearization_, 'newton') ...
                || strcmp(weno_linearization_, 'newton-FD-approx'), ...
                'weno_linearization option not possible. See available options');

            obj.weno_linearization = weno_linearization_;
        end
        
        function set.average_linearized_interface_wave_speed(obj, average_linearized_interface_wave_speed_)
           obj.average_linearized_interface_wave_speed = average_linearized_interface_wave_speed_;
        end
        
        function set.fd_weno_epsilon(obj, fd_weno_epsilon_)
           obj.fd_weno_epsilon = fd_weno_epsilon_;
        end
    
    end
    
end