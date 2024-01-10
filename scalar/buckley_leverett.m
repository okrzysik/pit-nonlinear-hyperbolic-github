% Class implementing the Buckley--Leverett equation
classdef buckley_leverett < cons_law_scalar 

    properties

    end
    
    methods

        % Constructor
        function obj = buckley_leverett(pde_pa_)
            % Call parent's constructor
            obj = obj@cons_law_scalar(pde_pa_, false);
           
        end
        
        function f = flux(obj, u)
           f = 4*u.^2 ./ (4*u.^2 + (1-u).^2);
        end
        
        function fprime = flux_jacobian(obj, u)
           fprime = 8*u.*(1-u) ./ ( 4*u.^2 + (1-u).^2 ).^2;
        end
        
        function id = id(obj)
            id = 'buckley-leverett';
        end
        
        function id = id_abbreviation(obj)
            id = 'BL';
        end
        
        
        % Evaluate the initial condition. 
        function u0 = initial_condition(obj, x)

            switch obj.pde_pa.ic_id

                case 3
                    u0 = zeros(size(x));
                    I = find(x > -1/2 & x < 0);
                    u0(I) = 1;

                case 6
                    % NOTE: This produces a shock for BL but no rarefaction.
                    u0 = 0.8 + 0.2*cos(pi*x); 
                    
                otherwise
                    error('ic_id = %d not recognised', obj.ic_id)
                     
            end
        
        % End of initial condition function
        end
        
        % Compute the maximum value of |fprime(z)| for z \in (min(u,v), max(u,v))
        % The optimization problem is complicated to solve in general because
        % the flux is not convex. Here we make the assumption that the
        % solution is in the interval [0,1], as is physically consistent
        % with the model of the PDE, and this simplifies the optimization
        % problem because there are only 3 choices for the maximum.
        function lambda = nonconvex_max_wave_speed(obj, u, v)
            
            %norm(u, inf)
            maxtol = 1e-6;
            %maxtol = 1;
            assert(norm(u, inf) <= 1+maxtol, 'sol violate max: 1-sol = %.2e\n', abs(1-norm(u, inf)))
            
            mintol = 1e-6;
            %mintol = 0;
            assert(min(u) >= 0-mintol, 'sol violate min: sol = %.2e\n', min(u))

            global_max     = 2.332030375854269;
            global_max_loc = 0.287140725416741;
            
            u_sorted = sort([u, v], 2); % Sort along each row, this is what the 2 does.
            u_left   = u_sorted(:, 1);
            u_right  = u_sorted(:, 2);
            
            % Initialize solution array.
            lambda = zeros(size(u));
            
            % Find pairs that have the location of the global max in
            % between them. For these pairs, the max of |fprime| over the
            % associated interval is just the global maximum.
            I = find(u_left < global_max_loc & global_max_loc < u_right);
            lambda(I) = global_max;
            
            J = setdiff(1:size(u_left), I);
            % Pairs for which both points are on one side or 
            % the other of the global maximum, the function is convex over
            % the associated interval, so we can just do a max.
            lambda(J) = max(obj.flux_jacobian(u_left(J)), obj.flux_jacobian(u_right(J)));
            
        end
        
    % End of methods.    
    end

% End of class.
end