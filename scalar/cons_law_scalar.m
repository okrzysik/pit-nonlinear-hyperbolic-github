% Abstract class implementing the 1D conservation law u_t + f(u)_x = 0.
%
% Uses Lax--Friedrichs numerical flux with the FV method.
% If a 1st-order FV method is used, then ERK1 is used in time, otherwise
% ERK3 is used in time. 
%
%
classdef cons_law_scalar < handle

    properties
        % Parameter structs
        pde_pa
        mesh_pa = []; 
        disc_pa = [];
        num_RK_stages = [];
        
        convex % Flag whether flux is convex.
        
        % Extremal values of problem over the spatial domain at time zero
        f0_prime_max = []; 
        u0_max       = [];
        u0_min       = [];
        
        linearization_pa   = []; % Parameters relating to the linearization
        linearization_data = []; % Data for linearizations
        
        % Class used to perform reconstructions
        reconstruction     = [];
    end
    
    methods

        % Constructor
        function obj = cons_law_scalar(pde_pa_, convex_)
            obj.pde_pa = pde_pa_;
            obj.convex = convex_;
        end

        function set.mesh_pa(obj, mesh_pa_)
            obj.mesh_pa  = mesh_pa_;
        end
        
        function set.disc_pa(obj, disc_pa_)
            obj.disc_pa  = disc_pa_;
            
            if disc_pa_.spatial_order == 1
                obj.num_RK_stages = 1;
            else 
                obj.num_RK_stages = 3;
            end
        end
        
        function set.reconstruction(obj, reconstruction_)
            assert(isa(reconstruction_, 'weighted_reconstruction'), ...
                'reconstruction required to be instance of weighted_reconstruction')
            
            obj.reconstruction = reconstruction_;
        end
        
        % Compute maximum |wave-speed| of equation over space-time based on the PDE
        % solution obeying a maximum principle and the range of the initial
        % conditon. 
        function [f0_prime_max, u0_min, u0_max] = compute_and_set_zero_extremal_values(obj, xmin, xmax)
            [u0_min, u0_max] = obj.extremal_u0(xmin, xmax);
            f0_prime_max = obj.max_abs_wave_speed(u0_min, u0_max, ~true);
            
            obj.u0_min       = u0_min;
            obj.u0_max       = u0_max;
            obj.f0_prime_max = f0_prime_max;
        end

        function set.linearization_data(obj, linearization_data_)
            assert(isa(linearization_data_, 'scalar_linearization_data'), ...
                'linearization_data required to be instance of linearization_data')

            obj.linearization_data = linearization_data_;
        end
        
        function set.linearization_pa(obj, linearization_pa_)
            assert(isa(linearization_pa_, 'scalar_linearization_pa'), ...
                'linearization_pa required to be instance of scalar_linearization_pa')
            
            obj.linearization_pa = linearization_pa_;
        end
        

        % Compute the maximum of |fprime(u)| on the interval [u_min, u_max].
        % Be careful here, it's not enough to comptute the discretize maxmimum of fprime on
        % the discretized initial condition because it can be zero on those
        % (e.g., as in the case for the BL flux when u0 contains only 0s 
        % and 1s.)
        function abs_fprime_max = max_abs_wave_speed(obj, u_min, u_max, plot_fluxes)
            
            [u0_maximizer, abs_fprime_min, EXITFLAG] = fminbnd(@(u) -abs(obj.flux_jacobian(u)), u_min, u_max);
            if EXITFLAG ~= 1
                error('fminbnd did not do what you intended');
            end
            abs_fprime_max = -abs_fprime_min;
            
            if plot_fluxes
                u0_range = linspace(u_min, u_max);
                figure()
                plot(u0_range, obj.flux(u0_range), '--m', 'LineWidth', 2, 'DisplayName', '$f(u)$')
                hold on
                plot(u0_range, fprime(u0_range), '-b', 'LineWidth', 2, 'DisplayName', '$f''(u)$')
                plot(u0_maximizer, abs_fprime_max, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 12, 'DisplayName', '$\max|f''(u)|$')
                hold off
                xlabel('$u$')
                lh = legend();
                lh.set('Location', 'Best') 
                title('maximum $|$wave-speed$|$ calculation')
                pause
            end
        end
        
        % Compute minimum and maximum of initial condition u0 as a function
        % of continuous x.
        function [u0_min, u0_max] = extremal_u0(obj, xmin, xmax)
            [~, u0_min, EXITFLAG1] = fminbnd(@(x)  obj.initial_condition(x), xmin, xmax);
            [~, u0_max, EXITFLAG2] = fminbnd(@(x) -obj.initial_condition(x), xmin, xmax);
            u0_max = -u0_max;
            if EXITFLAG1 ~= 1 || EXITFLAG2 ~= 1
                error('fminbnd did not do what you intended');
            end
        end
        
        
        
        %% step
        % Evolve u0 into u1
        function u1 = step(obj, t0_idx, u0)
            if ~isempty(obj.linearization_data); obj.linearization_data.t0_idx = t0_idx; end
            
            if obj.disc_pa.spatial_order == 1
                if ~isempty(obj.linearization_data); obj.linearization_data.stage_idx = 1; end
                u1 = u0 + obj.mesh_pa.dt * obj.spatial_discretization(u0);

            else
                if ~isempty(obj.linearization_data); obj.linearization_data.stage_idx = 1; end
                u0_1 =                  ( u0   + obj.mesh_pa.dt * obj.spatial_discretization(u0) );
                
                if ~isempty(obj.linearization_data); obj.linearization_data.stage_idx = 2; end
                u0_2 = 3/4 * u0 + 1/4 * ( u0_1 + obj.mesh_pa.dt * obj.spatial_discretization(u0_1) );
                
                if ~isempty(obj.linearization_data); obj.linearization_data.stage_idx = 3; end
                u1   = 1/3 * u0 + 2/3 * ( u0_2 + obj.mesh_pa.dt * obj.spatial_discretization(u0_2) ); 
            end

            if norm(u1, inf) > 10*norm(u0, inf)
                error('step: u1 > 10*u0')
            end
        end
        
        %% Linearized step 
        % Evolve e0 into e1
        function e1 = step_linearized(obj, t0_idx, e0)

            assert(~isempty(obj.linearization_data), 'Taking a linearized step requires non-empty linearization_data');
            obj.linearization_data.t0_idx = t0_idx;
            
            if obj.disc_pa.spatial_order == 1
                obj.linearization_data.stage_idx = 1;
                e1 = e0 + obj.mesh_pa.dt * obj.spatial_discretization_linearized(e0);

            else
                obj.linearization_data.stage_idx = 1;
                e0_1 =                  ( e0   + obj.mesh_pa.dt * obj.spatial_discretization_linearized(e0) );
                
                obj.linearization_data.stage_idx = 2;
                e0_2 = 3/4 * e0 + 1/4 * ( e0_1 + obj.mesh_pa.dt * obj.spatial_discretization_linearized(e0_1) );
                
                obj.linearization_data.stage_idx = 3;
                e1   = 1/3 * e0 + 2/3 * ( e0_2 + obj.mesh_pa.dt * obj.spatial_discretization_linearized(e0_2) ); 
            end

            if norm(e1, inf) > 10*norm(e0, inf)
                error('step_linearized: e1 > 10*e0')
            end
        end

        
        %% Spatial discretization
        % Evaluate the spatial discretization on current cell-averaged vector u0_bar.
        function L = spatial_discretization(obj, u0_bar)

            % Check if data needs to be stored for linearization
            linearize = ~isempty(obj.linearization_data);
            
            
            %% Reconstruct u0 at interfaces
            if ~linearize
                [u_left, u_right] = obj.reconstruction.interface_reconstruction(u0_bar);
            % The linearization procedure requires the reconstruction
            % weights to be stored.
            else
                [u_left, u_right, reconstruction_weights] = obj.reconstruction.interface_reconstruction(u0_bar);
            end

            %% Limit reconstructions if required
            if obj.disc_pa.limit_reconstructions
                u_left(u_left   < obj.u0_min) = obj.u0_min;
                u_right(u_right < obj.u0_min) = obj.u0_min;

                u_left(u_left   > obj.u0_max) = obj.u0_max;
                u_right(u_right > obj.u0_max) = obj.u0_max;
            end
            
            

            nx = obj.mesh_pa.nx;

        %     % Evaluation of the numerical flux.    
        %     f_num = pa.num_flux_nonlin;
        %     
        %     f_hat = [f_num(u_right(nx),     u_left(1)); ...    % First interface in the domain
        %              f_num(u_right(1:nx-1), u_left(2:nx)); ... % All interior interfaces in the domain
        %              f_num(u_right(nx),     u_left(1))];       % Last interface in the domain


            %% Evaluation of numerical flux on reconstructed u0
            % u_left and u_right are nx-dimensional arrays, holding the left and
            % right-hand side reconstructions for all nx cells. Convert them into
            % interface based arrays, so that u_neg and u_pos are nx+1-dimensional
            % arrays hold for every interface the reconstruction from below and the
            % reconstruction from above.

            % Reconstructions from below an interface, i.e., those that do
            % reconstruction in the right-hand side of a cell.
            % The reconstruction from below in the first interface is the 
            % reconstruction at the right interface of the final
            % cell.
            u_neg = [u_right(end); u_right];

            % Reconstructions from above an interface, i.e., those that do
            % reconstruction in the left-hand side of a cell.
            % The reconstruction from above in the last interface is the 
            % reconstruction done in the left-hand side of the first cell.
            u_pos  = [u_left; u_left(1)];

            
            % Determine the dissipation coefficient(s) nu.
            if strcmp(obj.disc_pa.num_flux_id, 'GLF')
                nu = obj.f0_prime_max;

            elseif strcmp(obj.disc_pa.num_flux_id, 'LLF')

                % nu is just the max of fprime evaluated at the two reconstructions
                % on each interface.
                if obj.convex
                    % Eval. wave-speed at interfaces. 
                    f_prime = struct('neg', obj.flux_jacobian(u_neg), 'pos', obj.flux_jacobian(u_pos));
                    
                    nu = max( abs(f_prime.neg), abs(f_prime.pos) );

                % Call a user-implemented max wave-speed function.
                elseif ~obj.convex
                    nu = obj.nonconvex_max_wave_speed(u_neg, u_pos);
                end

            end

            f_hat = 0.5*( (obj.flux(u_pos) + obj.flux(u_neg) ) - nu.*(u_pos - u_neg) );

            % f_hat_delta(i) = f_hat(i+1/2) - f_hat(i-1/2).
            f_hat_delta = f_hat(2:nx+1)-f_hat(1:nx);

            L = -1 / obj.mesh_pa.h * f_hat_delta;
            
            
            % Store data needed for linearization
            if linearize
                
                t0_idx    = obj.linearization_data.t0_idx;
                stage_idx = obj.linearization_data.stage_idx;
                
                % Eval wave-speed at interfaces if it wasn't already
                if ~obj.convex || strcmp(obj.disc_pa.num_flux_id, 'GLF')
                    f_prime = struct('neg', obj.flux_jacobian(u_neg), 'pos', obj.flux_jacobian(u_pos));
                end
                
                % Store the data
                obj.linearization_data.cell_averages{t0_idx}{stage_idx}          = u0_bar;
                obj.linearization_data.reconstruction_weights{t0_idx}{stage_idx} = reconstruction_weights;
                obj.linearization_data.dissipation{t0_idx}{stage_idx}            = nu;
                obj.linearization_data.wave_speed{t0_idx}{stage_idx}             = f_prime;
            end

        end
        % End of spatial discretization

        
        %% Linearized spatial discretization
        % Evaluate the spatial discretization on current cell-averaged vector e0.
        function L = spatial_discretization_linearized(obj, e0_bar)

            % Unpack parameters that tell us which linearization data to use.
            t0_idx    = obj.linearization_data.t0_idx;
            stage_idx = obj.linearization_data.stage_idx;
            
            %% Reconstruction of e0 at interfaces
            % How was the nonlinear reconstruction done? If linear then we use the
            % same linear reconstruction here.
            % Otherwise if a WENO reconstruction was used, then there are various
            % linearization options that can be used here. 
            if strcmp(obj.reconstruction.id, 'linear')
                [e0_left, e0_right] = obj.reconstruction.interface_reconstruction(e0_bar);

            elseif strcmp(obj.reconstruction.id, 'WENO')

                % A reconstruction with the weights used in the nonlinear reconstruction.
                if strcmp(obj.linearization_pa.weno_linearization, 'picard')
                    [e0_left, e0_right] = obj.reconstruction.interface_reconstruction(...
                        e0_bar, ...
                        obj.linearization_data.reconstruction_weights{t0_idx}{stage_idx}.left, ...
                        obj.linearization_data.reconstruction_weights{t0_idx}{stage_idx}.right);

                % Use exact gradient of WENO weights.    
                elseif strcmp(obj.linearization_pa.weno_linearization, 'newton')         
                    [e0_left, e0_right] = obj.reconstruction.interface_reconstruction_WENO_linearized(...
                        obj.linearization_data.cell_averages{t0_idx}{stage_idx}, ...
                        e0_bar);

                % Use a finite-difference approximation to evaluate the linearized
                % reconstruction.
                elseif strcmp(obj.linearization_pa.weno_linearization, 'newton-FD-approx')
                    [u0_left, u0_right] = obj.reconstruction.interface_reconstruction(...
                        obj.linearization_data.cell_averages{t0_idx}{stage_idx});

                    % It seems like the roundoff error can really mess things up
                    % here... I guess there is a linear and a nonlinear component
                    % when we expand out by the product rule. I think for any value
                    % of epsilon I get exactly the gradient of the linear term, and
                    % I think this one is most important as we approach
                    % convergence. 
                    [u0_eps_left, u0_eps_right] = obj.reconstruction.interface_reconstruction(...
                        obj.linearization_data.cell_averages{t0_idx}{stage_idx} + obj.linearization_pa.fd_weno_epsilon*e0_bar);

                    e0_left  = (u0_eps_left  - u0_left ) / obj.linearization_pa.fd_weno_epsilon;
                    e0_right = (u0_eps_right - u0_right) / obj.linearization_pa.fd_weno_epsilon;

                end
            end


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

            % Unpack linearized dissipation and wave-speeds
            nu        = obj.linearization_data.dissipation{t0_idx}{stage_idx};
            alpha_neg = obj.linearization_data.wave_speed{t0_idx}{stage_idx}.neg;
            alpha_pos = obj.linearization_data.wave_speed{t0_idx}{stage_idx}.pos;

            % Evaluate numerical flux
            f_hat = 0.5*( alpha_neg.*e0_neg + alpha_pos.*e0_pos + nu.*(e0_neg - e0_pos) );

            % f_hat_delta(i) = f_hat(i+1/2) - f_hat(i-1/2).
            f_hat_delta = f_hat(2:nx+1)-f_hat(1:nx);

            L = -1 / obj.mesh_pa.h * f_hat_delta;

        end
        % End of spatial_discretization_linearized
        
         
        
        
        
        %% Wrapper functions that must be implemented by child class
        function u0 = initial_condition(obj, x)
           error('Child class must implement "initial_condition"') 
        end

        function f = flux(obj, u)
           error('Child class must implement "flux"') 
        end
        
        function fprime = flux_jacobian(obj, u)
           error('Child class must implement "flux_jacobian"') 
        end
        
        % The local Lax--Friedrichs flux requires this for non-convex
        % fluxes. In principle it is possible to use the max wave speed
        % function implemented in this class. However, this is really slow
        % in general, so the user should really implement their own.
        % This function must compute the |largest| wave-speed between the
        % left and right states u and v, inclusive. 
        function lambda = nonconvex_max_wave_speed(obj, u, v)
            error('Child class must implement "nonconvex_max_wave_speed"') 
        end
        
        function id = id(obj)
            error('Child class must implement "PDE_id"') 
        end
        
        function id = id_abbreviation(obj)
            error('Child class must implement "PDE_id_abbreviation"') 
        end
        
    % End of methods.    
    end

% End of class.
end