% Abstract class implementing the 1D linear conservation law e_t +
% (alpha*e)_x = 0 where alpha is some prescribed function of x and t.
%
% A bunch of wave-speeds alpha are implemented here, with the specific one
% chosen by parsing the "wave_speed_id" to the constructor of this class.
%
% Note that this child class re-implements the step function from the
% parent class since that is not really equipped to deal with the linear 
% flux f(u) = alpha*u. Certain other functions are re-implemented also.
%
% Uses Lax--Friedrichs numerical flux with the FV method.
% If a 1st-order FV method is used, then ERK1 is used in time, otherwise
% ERK3 is used in time. 
%
classdef cons_lin_scalar < cons_law_scalar 

    properties
        wave_speed_id  = [];
        
        wave_speed_handle = []; % This is a handle to alpha as a function of x at some given time.
    end
    
    methods

        % Constructor
        function obj = cons_lin_scalar(pde_pa_, wave_speed_id_)
            % Call parent's constructor
            obj = obj@cons_law_scalar(pde_pa_, []);
            obj.wave_speed_id = wave_speed_id_;
        end
        
        % Wave-speeds for linear advection problems
        function a = wave_speed(obj, x, t)
            
            if obj.wave_speed_id == 1
                a = ones(size(x)); 
            elseif obj.wave_speed_id == 2
                a = -1/pi*sin(pi*x);
            elseif obj.wave_speed_id == 3
                a = 1/2*(1 + cos(pi*x).^2);
            elseif obj.wave_speed_id == 4
                a = -(sin(pi*(x-t)).^2);
            elseif obj.wave_speed_id == 5
                a = cos(2*pi*x) * cos(2*pi*t);
            else
               error('wave_speed_id not recognised') 
            end
            
        end
        
        % Re-implemented from parent class.
        function [f0_prime_max, u0_min, u0_max] = compute_and_set_zero_extremal_values(obj, xmin, xmax)
            obj.u0_min       = [];
            obj.u0_max       = [];
            obj.f0_prime_max = obj.max_abs_wave_speed();
        end

        % Re-implement from parent class. Max of abs wave-speed over
        % space-time.
        function abs_fprime_max = max_abs_wave_speed(obj)
           
            abs_fprime_max = 1;
            
            if obj.wave_speed_id == 2
                abs_fprime_max = 1 / pi;
            end
           
        end
                
        
        %% step
        % Evolve u0 into u1
        function e1 = step(obj, t0_idx, e0)
            
            if obj.disc_pa.spatial_order == 1
                t0 = obj.mesh_pa.t(t0_idx);
                obj.wave_speed_handle = @(x) obj.wave_speed(x, t0);
                e1 = e0 + obj.mesh_pa.dt * obj.spatial_discretization(e0);

            else
                
                c0 = 0;
                c1 = 1;
                c2 = 1/2;
                
                t0 = obj.mesh_pa.t(t0_idx);
                
                obj.wave_speed_handle = @(x) obj.wave_speed(x, t0 + c0 * obj.mesh_pa.dt);
                e0_1 =                  ( e0   + obj.mesh_pa.dt * obj.spatial_discretization(e0) );
                
                obj.wave_speed_handle = @(x) obj.wave_speed(x, t0 + c1 * obj.mesh_pa.dt);
                e0_2 = 3/4 * e0 + 1/4 * ( e0_1 + obj.mesh_pa.dt * obj.spatial_discretization(e0_1) );
                
                obj.wave_speed_handle = @(x) obj.wave_speed(x, t0 + c2 * obj.mesh_pa.dt);
                e1   = 1/3 * e0 + 2/3 * ( e0_2 + obj.mesh_pa.dt * obj.spatial_discretization(e0_2) ); 
            end

            if norm(e1, inf) > 10*norm(e0, inf)
                error('step: u1 > 10*u0')
            end
        end
        
       
        %% Spatial discretization
        % Evaluate the spatial discretization on current cell-averaged vector e0.
        % This uses a LF flux where the interfacial wave-speeds are
        % evaluated directly rather than being reconstructed.
        function L = spatial_discretization(obj, e0_bar)
            
            %% Reconstruction of e0 at interfaces
            [e0_left, e0_right] = obj.reconstruction.interface_reconstruction(e0_bar);


            %% Evaluation of the numerical flux on reconstructed e0 
            nx = obj.mesh_pa.nx;

            % e_left and e_right are nx-dimensional arrays, holding the left and
            % right-hand side reconstructions for all nx cells. Convert them into
            % interface based arrays, so that u_neg and u_pos are nx+1-dimensional
            % arrays hold for every interface the reconstruction from below and the
            % reconstruction from above.

            % Reconstructions from below an interface, i.e., those that do
            % reconstruction in the right-hand side of a cell.
            % The reconstruction from below in the first interface is the 
            % reconstruction at the right interface of the final
            % cell.
            e0_neg = [e0_right(end); e0_right];

            % Reconstructions from above an interface, i.e., those that do
            % reconstruction in the left-hand side of a cell.
            % The reconstruction from above in the last interface is the 
            % reconstruction done in the left-hand side of the first cell.
            e0_pos  = [e0_left; e0_left(1)];

            % Evaluate the wave-speed function on all FV cell interfaces
            alpha_interfaces = obj.wave_speed_handle(obj.mesh_pa.x_interfaces);
            
            % For all interfaces, get wave-speed from above and below (we
            % don't really do this because we don't reconstruct the
            % wave-speed, we just evaluate it at each interface).
            alpha_neg = alpha_interfaces;
            alpha_pos = alpha_interfaces;

            % Determine the dissipation coefficient(s) nu.
            if strcmp(obj.disc_pa.num_flux_id, 'GLF')
                nu = obj.max_abs_wave_speed;

            elseif strcmp(obj.disc_pa.num_flux_id, 'LLF')
                nu = abs(alpha_interfaces);

            end
            
            % Evaluate numerical flux
            f_hat = 0.5*( alpha_neg.*e0_neg + alpha_pos.*e0_pos + nu.*(e0_neg - e0_pos) );

            % f_hat_delta(i) = f_hat(i+1/2) - f_hat(i-1/2).
            f_hat_delta = f_hat(2:nx+1)-f_hat(1:nx);

            L = -1 / obj.mesh_pa.h * f_hat_delta;

        end
        % End of spatial_discretization
           
        
        %% Wrapper functions that must be implemented by this class
        function u0 = initial_condition(obj, x)
           switch obj.pde_pa.ic_id
                % moving shock    
                case 1
                    u0 = sin(pi*x).^4;

                case 2    
                    u0 = -sin(pi*x); 

                case 3
                    u0 = zeros(size(x));
                    I = find(x > -1/2 & x < 0);
                    u0(I) = 1;
                    
                otherwise
                    error('ic_id = %d not recognised', obj.ic_id)         
            end 
        end

        function id = id(obj)
            id = 'linear';
        end
        
        function id = id_abbreviation(obj)
            id = 'L';
        end
        
    % End of methods.    
    end

% End of class.
end