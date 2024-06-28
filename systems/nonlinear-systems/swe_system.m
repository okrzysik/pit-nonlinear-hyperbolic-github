% Class implementing shallow water equations in 1D.
%
% Where spatial vectors are parsed into and out of functions they are
% expected to be blocked by variable.
%
% Numerical flux can be either Roe or local Lax--Friedrichs. Set in
% disc_pa.num_flux_id as either "ROE" or "LLF"
%
classdef swe_system < cons_law_system
   
    properties

    end
    
   
    methods
        
    %% Constructor
    function obj = swe_system(pde_pa_)  
        
        % Call parent's constructor
        obj = obj@cons_law_system(pde_pa_, 2);
        
    end
    
    % Evaluate the spatial discretization on current cell-averaged vector q.
    % There's the option to output some linearization data.
    function L = spatial_discretization(obj, t0idx, q0)

        nx = obj.mesh_pa.nx;
        % Extract q1 and q2 from stacked q vector.
        q1 = q0(1:nx);
        assert(min(q1) >= 0, 'SWE height is negative at t0idx = %d', t0idx)
        q2 = q0(nx+1:2*nx);

        % Reconstructions at left and right interfaces of all nx cells. 
        % Reconstructions are 1st-order: The solution is assumed constant in 
        % each cell. More generally, we could use high-order WENO
        % reconstructions here.
        q1_left  = q1;
        q1_right = q1;
        q2_left  = q2;
        q2_right = q2;

        % Get values in single ghost cell at each end of the domain. 
        % Periodic boundaries: The right-hand side ghost cell is the same
        % as the first physical cell. The left-hand side ghost cell is the
        % same as the last physical cell
        if strcmp(obj.pde_pa.bcs, 'periodic')
            h_rbc  = q1_left(1);
            h_lbc  = q1_right(nx);
            hu_rbc = q2_left(1);
            hu_lbc = q2_right(nx);
            
        % Constant extrapolation boundary conditions. This means the
        % boundary is just set whatever is in the cell next to it.
        elseif strcmp(obj.pde_pa.bcs, 'constant')
            h_rbc  = q1_right(end);
            h_lbc  = q1_left(1);
            hu_rbc = q2_right(end);
            hu_lbc = q2_left(1);
            
        end

        % Re-order and extend arrays to give reconstructions at all nx+1
        % interfaces, this uses periodicity since the nx+1st interface is the
        % same as the 1st interface.
        % "neg" means from below an interface, i.e., the reconstruction in the 
        % RHS of the cell to the interface's left.
        % "pos" means from above an interface, i.e., the reconstruction in the 
        % LHS of the cell to the interface's right.
        q1_pos = [q1_left; h_rbc   ];
        q1_neg = [h_lbc;   q1_right];
        q2_pos = [q2_left; hu_rbc  ];
        q2_neg = [hu_lbc;  q2_right];

        % Compute primitive variables
        % q1 = h, so h = q1.
        h_pos = q1_pos; 
        h_neg = q1_neg;
        % q2 = h*u, so u = q2 / h.
        u_pos = q2_pos ./ h_pos; 
        u_neg = q2_neg ./ h_neg;
        
        if ~isempty(obj.linearization_data)
            obj.linearization_data.h_neg{obj.linearization_data.t0_idx} = h_neg;
            obj.linearization_data.h_pos{obj.linearization_data.t0_idx} = h_pos;
            obj.linearization_data.u_neg{obj.linearization_data.t0_idx} = u_neg;
            obj.linearization_data.u_pos{obj.linearization_data.t0_idx} = u_pos;
        end

        % Compute physical flux at all interfaces
        [f1_pos, f2_pos] = obj.flux(h_pos, u_pos);
        [f1_neg, f2_neg] = obj.flux(h_neg, u_neg);

        %% Roe flux
        if strcmp(obj.disc_pa.num_flux_id, 'ROE')

            % Compute Roe averages all all interfaces. See LeVeque p. 321.
            h_hat = 0.5*(h_neg + h_pos);
            u_hat = (sqrt(h_neg).*u_neg + sqrt(h_pos).*u_pos) ./ (sqrt(h_neg) + sqrt(h_pos));
            c_hat = sqrt(h_hat);

            % The jump for each conserved variable: 
            delta1 = q1_pos - q1_neg;
            delta2 = q2_pos - q2_neg;

            % Coefficients representing jump in eigenvector basis.
            alpha1 = (  (u_hat + c_hat) .* delta1 - delta2 ) ./ (2*c_hat);
            alpha2 = ( -(u_hat - c_hat) .* delta1 + delta2 ) ./ (2*c_hat);

            % Absolute eigenvalues
            %lambda1_abs = abs(u_hat - c_hat);
            %lambda2_abs = abs(u_hat + c_hat);
            % Smoothed absolute value, as in Harten's entropy fix (see
            % LeVeque p. 326).
            lambda1_abs = obj.smoothed_absolute_value(u_hat - c_hat);
            lambda2_abs = obj.smoothed_absolute_value(u_hat + c_hat);

            % Numerical flux
            % Note the eigenvector components for q1 are just 1, so no need to
            % multiply by them.
            f_hat1 = 0.5*( f1_pos + f1_neg - ...
                (lambda1_abs.*alpha1 + lambda2_abs.*alpha2 ));
            % Eigenvector components for q2 are u+c and u-c
            f_hat2 = 0.5*( f2_pos + f2_neg - ...
                (lambda1_abs.*alpha1.*(u_hat - c_hat) + lambda2_abs.*alpha2.*(u_hat + c_hat) ));

        %% Lax--Friedrichs
        elseif strcmp(obj.disc_pa.num_flux_id, 'LLF')

            % Take max of max of 2 eigenvalues. 
            % (I assume this requires the flux to be "convex" as in the
            % scalar case?)
            lambda1_max = max( abs(u_pos - sqrt(h_pos)), abs(u_neg - sqrt(h_neg)) ); % u - sqrt(h)
            lambda2_max = max( abs(u_pos + sqrt(h_pos)), abs(u_neg + sqrt(h_neg)) ); % u + sqrt(h)
            lambda_max  = max(lambda1_max, lambda2_max);

            % The jump for each conserved variable: 
            delta1 = q1_pos - q1_neg;
            delta2 = q2_pos - q2_neg;

            % Numerical flux
            f_hat1 = 0.5*( f1_pos + f1_neg - lambda_max .* delta1 );
            f_hat2 = 0.5*( f2_pos + f2_neg - lambda_max .* delta2 );
        end

        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat1(2:nx+1) - f_hat1(1:nx); ...     % First component.
               f_hat2(2:nx+1) - f_hat2(1:nx)] / obj.mesh_pa.h; % Second component.
    end
    % End of shallow_water_spatial_discretization
    
    
    
    % Evaluate linearized spatial discretization on current cell-averaged vector e.
    function L = spatial_discretization_linearized(obj, t0_idx, e0)

        %assert(~strcmp(obj.disc_pa.num_flux_id, 'LLF'), 'Linearized LLF not implemented')
        
        nx = obj.mesh_pa.nx;
        %% Compute reconstructions of linear variables
        % Extract e1 and e2 from stacked e vector.
        e1 = e0(1:nx);
        e2 = e0(nx+1:2*nx);

        % Reconstructions at left and right interfaces of all nx cells. 
        % Reconstructions are 1st-order: The solution is assumed constant in 
        % each cell. More generally, we could use high-order WENO
        % reconstructions here.
        e1_left  = e1;
        e1_right = e1;
        e2_left  = e2;
        e2_right = e2;

        % Get values in single ghost cell at each end of the domain. 
        % Periodic boundaries: The right-hand side ghost cell is the same
        % as the first physical cell. The left-hand side ghost cell is the
        % same as the last physical cell
        if strcmp(obj.pde_pa.bcs, 'periodic')
            e1_rbc = e1_left(1);
            e1_lbc = e1_right(nx);
            e2_rbc = e2_left(1);
            e2_lbc = e2_right(nx);
            
        % Constant extrapolation boundary conditions for nonlinear problem,
        % I believe this means zero BCs for the error... But not totally
        % sure...
        elseif strcmp(obj.pde_pa.bcs, 'constant')
            e1_rbc = 0;
            e1_lbc = 0;
            e2_rbc = 0;
            e2_lbc = 0;

            e1_rbc = e1(end);
            e1_lbc = e1(1);
            e2_rbc = e2(end);
            e2_lbc = e2(1);
            
        end

        e1_pos = [e1_left; e1_rbc];
        e1_neg = [e1_lbc;  e1_right];
        e2_pos = [e2_left; e2_rbc];
        e2_neg = [e2_lbc;  e2_right];
        
        h_pos = obj.linearization_data.h_pos{obj.linearization_data.t0_idx};
        h_neg = obj.linearization_data.h_neg{obj.linearization_data.t0_idx};
        u_pos = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};
        
        %% Compute the numerical flux
        % Compute physical flux at all interfaces
        [f1_pos, f2_pos] = obj.flux_linearized(e1_pos, e2_pos, h_pos, u_pos);
        [f1_neg, f2_neg] = obj.flux_linearized(e1_neg, e2_neg, h_neg, u_neg);

        % The jump for each conserved variable: 
        delta1 = e1_pos - e1_neg;
        delta2 = e2_pos - e2_neg;

        %% Roe flux
        if strcmp(obj.disc_pa.num_flux_id, 'ROE')

            h_hat = 0.5*(h_neg + h_pos);
            u_hat = (sqrt(h_neg).* u_neg + sqrt(h_pos).* u_pos) ./ (sqrt(h_neg) + sqrt(h_pos));
            c_hat = sqrt(h_hat);

            % Coefficients representing jump in eigenvector basis.
            beta1 = (  (u_hat + c_hat) .* delta1 - delta2 ) ./ (2*c_hat);
            beta2 = ( -(u_hat - c_hat) .* delta1 + delta2 ) ./ (2*c_hat);

            % Absolute eigenvalues
            %lambda1_abs = abs(u_hat - c_hat);
            %lambda2_abs = abs(u_hat + c_hat);
            % Smoothed absolute value, as in Harten's entropy fix (see
            % LeVeque p. 326).
            lambda1_abs = obj.smoothed_absolute_value(u_hat - c_hat);
            lambda2_abs = obj.smoothed_absolute_value(u_hat + c_hat);

            % Numerical flux
            % Note the eigenvector components for e1 are just 1, so no need to
            % multiply by them.
            f_hat1 = 0.5*( f1_pos + f1_neg - ...
                (lambda1_abs.*beta1 + lambda2_abs.*beta2 ));
            % Eigenvector components for q2 are u+c and u-c
            f_hat2 = 0.5*( f2_pos + f2_neg - ...
                (lambda1_abs.*beta1.*(u_hat - c_hat) + lambda2_abs.*beta2.*(u_hat + c_hat) ));

        %% Lax--Friedrichs
        elseif strcmp(obj.disc_pa.num_flux_id, 'LLF')

            % Take max of max of 2 eigenvalues. 
            % (I assume this requires the flux to be "convex" as in the
            % scalar case?)
            lambda1_max = max( abs(u_pos - sqrt(h_pos)), abs(u_neg - sqrt(h_neg)) ); % u - sqrt(h)
            lambda2_max = max( abs(u_pos + sqrt(h_pos)), abs(u_neg + sqrt(h_neg)) ); % u + sqrt(h)
            lambda_max  = max(lambda1_max, lambda2_max);

            % The jump for each conserved variable: 
            delta1 = e1_pos - e1_neg;
            delta2 = e2_pos - e2_neg;

            % Numerical flux
            f_hat1 = 0.5*( f1_pos + f1_neg - lambda_max .* delta1 );
            f_hat2 = 0.5*( f2_pos + f2_neg - lambda_max .* delta2 );
        end

        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat1(2:nx+1) - f_hat1(1:nx); ...     % First component.
               f_hat2(2:nx+1) - f_hat2(1:nx)] / obj.mesh_pa.h; % Second component.
    end
    % End of spatial_discretization
    
    function L = spatial_discretization_linearized_characteristic_variable(obj, t0_idx, w0, k)

        nx = obj.mesh_pa.nx;
        %% Compute reconstructions

        % Reconstructions at left and right interfaces of all nx cells. 
        % Reconstructions are 1st-order: The solution is assumed constant in 
        % each cell. More generally, we could use high-order WENO
        % reconstructions here.
        w_left  = w0;
        w_right = w0;

        % Get values in single ghost cell at each end of the domain. 
        % Periodic boundaries: The right-hand side ghost cell is the same
        % as the first physical cell. The left-hand side ghost cell is the
        % same as the last physical cell
        if strcmp(obj.pde_pa.bcs, 'periodic')
            w_rbc = w_left(1);
            w_lbc = w_right(nx);
            
        % Constant extrapolation boundary conditions for nonlinear problem,
        % I believe this means zero BCs for the error... But not totally
        % sure...
        elseif strcmp(obj.pde_pa.bcs, 'constant')
            w_rbc = 0;
            w_lbc = 0;

            w_rbc = w0(nx);
            w_lbc = w0(1);
            
        end

        w_pos = [w_left; w_rbc];
        w_neg = [w_lbc;  w_right];
        
        h_pos = obj.linearization_data.h_pos{obj.linearization_data.t0_idx};
        h_neg = obj.linearization_data.h_neg{obj.linearization_data.t0_idx};
        u_pos = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};
        
        % lambda1 = u - sqrt(h)
        if k == 1
            lambda_pos = u_pos - sqrt(h_pos);
            lambda_neg = u_neg - sqrt(h_neg);
        % lambda2 = u + sqrt(h)
        elseif k == 2
            lambda_pos = u_pos + sqrt(h_pos);
            lambda_neg = u_neg + sqrt(h_neg);
        else
            error('wave-speed index must be either 1 or 2.')
        end
        
        %% Compute the numerical flux
        % Compute physical flux at all interfaces
        f_pos = lambda_pos .* w_pos;
        f_neg = lambda_neg .* w_neg;

        % The jump for each conserved variable: 
        delta = w_pos - w_neg;

        
        %% Roe flux
        if strcmp(obj.disc_pa.num_flux_id, 'ROE')

            h_hat = 0.5*(h_neg + h_pos);
            u_hat = (sqrt(h_neg).* u_neg + sqrt(h_pos).* u_pos) ./ (sqrt(h_neg) + sqrt(h_pos));
            c_hat = sqrt(h_hat);

            % Absolute eigenvalues
            %lambda1_abs = abs(u_hat - c_hat);
            %lambda2_abs = abs(u_hat + c_hat);
            % Smoothed absolute value, as in Harten's entropy fix (see
            % LeVeque p. 326).
            if k == 1
                lambda_abs = obj.smoothed_absolute_value(u_hat - c_hat);
            elseif k == 2
                lambda_abs = obj.smoothed_absolute_value(u_hat + c_hat);
            end
            
            f_hat = 0.5*( f_pos + f_neg - lambda_abs .* delta );

        %% Lax--Friedrichs
        elseif strcmp(obj.disc_pa.num_flux_id, 'LLF')

           f_hat = 0.5*( f_pos + f_neg - max(abs(lambda_pos),abs(lambda_neg)) .* delta );
           %f_hat = 0.5*( f_pos + f_neg - 0.5*(abs(lambda_pos)+abs(lambda_neg)) .* delta );
        end

        
        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat(2:nx+1) - f_hat(1:nx)] / obj.mesh_pa.h; 
    end
    % End of spatial_discretization_linearized_characteristic_variable
    
    
    
    % SWE flux. See LeVeque p. 320
    function [f1, f2] = flux(obj, h, u)
        f1 = h .* u; 
        f2 = h .* u.^2 + 0.5*h.^2;
    end
    
    % Compute the linear flux f = A(q)*e, where A is the Jacobian of the
    % nonlinear flux f evaluated at q.
    % See LeVeque p. 255.
    function [f1, f2] = flux_linearized(obj, e1, e2, h, u)
        f1 = e2;
        f2 = (-u.^2 + h) .* e1 + 2 * u.*e2;
    end
    
    % ic_description returns the name of the IC along with any parameters
    function [h, u, ic_description] = initial_condition(obj, x)

        if strcmp(obj.pde_pa.ic_id, 'idp1')
            beta = 5;
            h = 1 + obj.pde_pa.ic_epsilon*exp(-beta * (x-2.5).^2);
            u = zeros(size(x));

            ic_description = sprintf('idp1-eps%.2f', obj.pde_pa.ic_epsilon);

        elseif strcmp(obj.pde_pa.ic_id, 'idp2')
            u = zeros(size(x));
            %h = 1 + 0.5*cos(pi*x/2);

            h = 1 + obj.pde_pa.ic_epsilon*cos(pi*x/5);
            
            ic_description = sprintf('idp2-eps%.2f', obj.pde_pa.ic_epsilon);
            
        elseif strcmp(obj.pde_pa.ic_id, 'dam-break')
            
            h = ones(size(x));
            if strcmp(obj.pde_pa.bcs, 'periodic')
                I = find(x > -2 & x < 2); % Periodic-bc-friendly version   
            else
                I = find(x < 0);
            end

            % h = 1 + eps, eps = 2 is the standard dam break problem from
            % LeVeque p. 259
            h(I) = 1 + obj.pde_pa.ic_epsilon; 
            u = zeros(size(x));

            ic_description = sprintf('dam-break-eps%.2f', obj.pde_pa.ic_epsilon);

        elseif strcmp(obj.pde_pa.ic_id, 'two-shock')
            h = ones(size(x));
            u = ones(size(x));

            I = find(x > 0);
            u(I) = -1;   
        end
    end
    % End of initial condition
    
    function lambda = wave_speeds(obj, q)
        nx = obj.mesh_pa.nx;
        h = q(1:nx);
        u = q(nx+1:2*nx) ./ h;
        
        lambda = [u - sqrt(h); u + sqrt(h)];
    end
    
    
    % Map from characteristic variables b to primitive variables a. This
    % is achieved by applying R to b.
    % q is the state used to evaluate the eigenvectors.
    function a = right_eigenvector_map(obj, q, b) 
        nx = obj.mesh_pa.nx;
        h = q(1:nx);
        u = q(nx+1:2*nx) ./ h;
        
        lambda1 = u - sqrt(h);
        lambda2 = u + sqrt(h);
        
        a            = zeros(size(b));
        a(1:nx)      =            b(1:nx) +            b(nx+1:2*nx);
        a(nx+1:2*nx) = lambda1 .* b(1:nx) + lambda2 .* b(nx+1:2*nx);
    end
    
    % Map from primitive variables a to characteristic variables b. This
    % is achieved by applying R^-1 to a.
    % q is the state used to evaluate the eigenvectors.
    function b = left_eigenvector_map(obj, q, a) 
        nx = obj.mesh_pa.nx;
        h = q(1:nx);
        u = q(nx+1:2*nx) ./ h;
        
        lambda1 = u - sqrt(h);
        lambda2 = u + sqrt(h);
        
        b            = zeros(size(a));
        b(1:nx)      = ( lambda2 .* a(1:nx) - a(nx+1:2*nx)) ./ (lambda2 - lambda1);
        b(nx+1:2*nx) = (-lambda1 .* a(1:nx) + a(nx+1:2*nx)) ./ (lambda2 - lambda1);
    end

    % Return the kth linearized wave-speed based at time t0 based on the
    % stored linearization data.
    % "neg" means from below an interface, i.e., the reconstruction in the 
    % RHS of the cell to the interface's left.
    % "pos" means from above an interface, i.e., the reconstruction in the 
    % LHS of the cell to the interface's right.
    function [lambda_k_neg, lambda_k_pos] = linearized_wave_speed(obj, k) 
        assert(~isempty(obj.linearization_data), 'Cannot compute linearized wave-speed without linearization data');

        h_pos = obj.linearization_data.h_pos{obj.linearization_data.t0_idx};
        h_neg = obj.linearization_data.h_neg{obj.linearization_data.t0_idx};
        u_pos = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};

        if k == 1
            lambda_k_neg = u_neg - sqrt( h_neg );
            lambda_k_pos = u_pos - sqrt( h_pos );
        elseif k == 2
            lambda_k_neg = u_neg + sqrt( h_neg );
            lambda_k_pos = u_pos + sqrt( h_pos );
        end
        
    end
    
    end
    % End of methods

end