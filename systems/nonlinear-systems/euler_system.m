% Class implementing Euler equations in 1D.
%
% Where spatial vectors are parsed into and out of functions they are
% expected to be blocked by variable.
%
% Numerical flux can be either Roe or local Lax--Friedrichs. Set in
% disc_pa.num_flux_id as either "ROE" or "LLF"
%
classdef euler_system < cons_law_system
   
    properties
        gamma = 7/5; % Adiabatic constant. Default value is 7/5. 
    end
    
   
    methods
        
    %% Constructor
    function obj = euler_system(pde_pa_)  
        
        % Call parent's constructor
        obj = obj@cons_law_system(pde_pa_, 3);
        
    end

    % Set value of gamma
    function obj = set.gamma(obj, gamma_)
        obj.gamma = gamma_;
    end
    
    % Evaluate the spatial discretization on current cell-averaged vector q.
    % There's the option to output some linearization data.
    function L = spatial_discretization(obj, t0_idx, q0)

        nx = obj.mesh_pa.nx;
    
        % Extract q1, q2, and q3 from stacked q vector.
        q1 = q0(1:nx);
        assert(min(q1) >= 0, 'Euler: Density is negative at t0idx = %d', t0_idx)
        assert(isreal(q1),   'Euler: Density is complex at t0idx = %d', t0_idx)
        q2 = q0(nx+1:2*nx);
        q3 = q0(2*nx+1:3*nx);

        % Reconstructions at left and right interfaces of all nx cells. 
        % Reconstructions are 1st-order: The solution is assumed constant in 
        % each cell. More generally, we could use high-order WENO
        % reconstructions here.
        q1_left  = q1;
        q1_right = q1;
        q2_left  = q2;
        q2_right = q2;
        q3_left  = q3;
        q3_right = q3;
        
        % Get values in single ghost cell at each end of the domain. 
        % Periodic boundaries: The right-hand side ghost cell is the same
        % as the first physical cell. The left-hand side ghost cell is the
        % same as the last physical cell
        if strcmp(obj.pde_pa.bcs, 'periodic')
            q1_rbc = q1_left(1);
            q1_lbc = q1_right(nx);
            q2_rbc = q2_left(1);
            q2_lbc = q2_right(nx);
            q3_rbc = q3_left(1);
            q3_lbc = q3_right(nx);
            
        % Constant extrapolation boundary conditions. This means the
        % boundary is just set whatever is in the cell next to it.
        elseif strcmp(obj.pde_pa.bcs, 'constant')
            q1_rbc = q1_right(end);
            q1_lbc = q1_left(1);
            q2_rbc = q2_right(end);
            q2_lbc = q2_left(1);
            q3_rbc = q3_right(end);
            q3_lbc = q3_left(1);
            
        end

        % Re-order and extend arrays to give reconstructions at all nx+1
        % interfaces, this uses periodicity since the nx+1st interface is the
        % same as the 1st interface.
        % "neg" means from below an interface, i.e., the reconstruction in the 
        % RHS of the cell to the interface's left.
        % "pos" means from above an interface, i.e., the reconstruction in the 
        % LHS of the cell to the interface's right.
        q1_pos = [q1_left; q1_rbc   ];
        q1_neg = [q1_lbc;  q1_right];
        q2_pos = [q2_left; q2_rbc  ];
        q2_neg = [q2_lbc;  q2_right];
        q3_pos = [q3_left; q3_rbc  ];
        q3_neg = [q3_lbc;  q3_right];
        

        % Compute primitive variables
        % q1 = rho, so rho = q1.
        rho_pos = q1_pos; 
        rho_neg = q1_neg;
        % q2 = rho*u, so u = q2 / rho.
        u_pos = q2_pos ./ rho_pos; 
        u_neg = q2_neg ./ rho_neg;
        % q3 = E, so E = q3
        E_pos = q3_pos; 
        E_neg = q3_neg;
        
        if ~isempty(obj.linearization_data)
            obj.linearization_data.rho_neg{obj.linearization_data.t0_idx} = rho_neg;
            obj.linearization_data.rho_pos{obj.linearization_data.t0_idx} = rho_pos;
            obj.linearization_data.u_neg{obj.linearization_data.t0_idx}   = u_neg;
            obj.linearization_data.u_pos{obj.linearization_data.t0_idx}   = u_pos;
            obj.linearization_data.E_neg{obj.linearization_data.t0_idx}   = E_neg;
            obj.linearization_data.E_pos{obj.linearization_data.t0_idx}   = E_pos;
        end

        % Compute physical flux at all interfaces
        [f1_pos, f2_pos, f3_pos, H_pos] = obj.flux(rho_pos, u_pos, E_pos);
        [f1_neg, f2_neg, f3_neg, H_neg] = obj.flux(rho_neg, u_neg, E_neg);

        % Compute Roe averages at all interfaces. See LeVeque p. 322--323.
        u_hat = (sqrt(rho_neg).*u_neg + sqrt(rho_pos).*u_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));
        H_hat = (sqrt(rho_neg).*H_neg + sqrt(rho_pos).*H_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));

        c_hat = sqrt( (obj.gamma-1)*( H_hat - 0.5*u_hat.^2 ) );

        %max(abs(c_hat + abs(u_hat)))

        % The jump for each conserved variable: 
        delta1 = q1_pos - q1_neg;
        delta2 = q2_pos - q2_neg;
        delta3 = q3_pos - q3_neg;

        % Coefficients representing jump in eigenvector basis.
        alpha2 = (obj.gamma - 1) .* ( (H_hat - u_hat.^2) .* delta1 + u_hat.*delta2 - delta3 ) ./ (c_hat.^2);

        % Based on LeVeque
        alpha3 = (delta2 + (c_hat - u_hat).*delta1 - c_hat.*alpha2) ./ (2*c_hat);
        alpha1 = delta1 - alpha2 - alpha3;

        % Based on Hesthaven. Are these the same? It would seem so.
    %     alpha11 = ( (c_hat + u_hat).*delta1 - delta2 - c_hat.*alpha2) ./ (2*c_hat);
    %     alpha33 = delta1 - alpha2 - alpha1;
    %     norm(alpha3 - alpha33)
    %     norm(alpha1 - alpha11)

        % Absolute eigenvalues
%         lambda1_abs = abs(u_hat - c_hat);
%         lambda2_abs = abs(u_hat);
%         lambda3_abs = abs(u_hat + c_hat);
        % Smoothed absolute value, as in Harten's entropy fix (see LeVeque p. 326).
        lambda1_abs = obj.smoothed_absolute_value(u_hat - c_hat);
        lambda2_abs = obj.smoothed_absolute_value(u_hat);
        lambda3_abs = obj.smoothed_absolute_value(u_hat + c_hat);

        % Numerical flux.
        % See LeVeque, p.300--301 for eigenvalues and vectors.
        % Note the eigenvector components for q1 are just 1, so no need to
        % multiply by them.
        f_hat1 = 0.5*( f1_pos + f1_neg - ...
            (lambda1_abs.*alpha1 + ...
             lambda2_abs.*alpha2 + ...
             lambda3_abs.*alpha3 ) );
        % Eigenvector components for q2 are u-c, u, and u+c
        f_hat2 = 0.5*( f2_pos + f2_neg - ...
            (lambda1_abs.*alpha1.*(u_hat - c_hat) + ...
             lambda2_abs.*alpha2.*(u_hat ) + ...
             lambda3_abs.*alpha3.*(u_hat + c_hat)) );
        % Eigenvector components for q3 are H-u*c, 1/2*u^2, and H + u*c
        f_hat3 = 0.5*( f3_pos + f3_neg - ...
            (lambda1_abs.*alpha1.*(H_hat - u_hat.*c_hat) + ...
             lambda2_abs.*alpha2.*(0.5*u_hat.^2 ) + ...
             lambda3_abs.*alpha3.*(H_hat + u_hat.*c_hat)) );

        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat1(2:nx+1) - f_hat1(1:nx); ...     % First component.
               f_hat2(2:nx+1) - f_hat2(1:nx); ...     % Second component.
               f_hat3(2:nx+1) - f_hat3(1:nx)] / obj.mesh_pa.h; % Third component.
    end
    % End of spatial_discretization
    
    % Euler flux in conservative variables.
    function [f1, f2, f3, H] = flux(obj, rho, u, E)
        p  = obj.pressure(rho, u, E);
        H  = (E + p) ./ rho;
        
        f1 = rho .* u; 
        f2 = rho .* u.^2 + p;
        f3 = (E + p) .* u;
    end
    
    % The ideal polytropic EOS. See LeVeque p. 295
    function p = pressure(obj, rho, u, E)
        p = (obj.gamma - 1) .* (E - 0.5*rho.*u.^2);
    end
    
    
    % Evaluate the spatial discretization on current cell-averaged vector q.
    % There's the option to output some linearization data.
    function L = spatial_discretization_linearized(obj, t0_idx, e0)

        nx = obj.mesh_pa.nx;
    
        % Extract q1, q2, and q3 from stacked q vector.
        e1 = e0(1:nx);
        e2 = e0(nx+1:2*nx);
        e3 = e0(2*nx+1:3*nx);

        % Reconstructions at left and right interfaces of all nx cells. 
        % Reconstructions are 1st-order: The solution is assumed constant in 
        % each cell. More generally, we could use high-order WENO
        % reconstructions here.
        e1_left  = e1;
        e1_right = e1;
        e2_left  = e2;
        e2_right = e2;
        e3_left  = e3;
        e3_right = e3;
        
        % Get values in single ghost cell at each end of the domain. 
        % Periodic boundaries: The right-hand side ghost cell is the same
        % as the first physical cell. The left-hand side ghost cell is the
        % same as the last physical cell
        if strcmp(obj.pde_pa.bcs, 'periodic')
            e1_rbc = e1_left(1);
            e1_lbc = e1_right(nx);
            e2_rbc = e2_left(1);
            e2_lbc = e2_right(nx);
            e3_rbc = e3_left(1);
            e3_lbc = e3_right(nx);
            
        % Constant extrapolation boundary conditions. This means the
        % boundary is just set whatever is in the cell next to it.
        elseif strcmp(obj.pde_pa.bcs, 'constant')
            e1_rbc = e1_right(end);
            e1_lbc = e1_left(1);
            e2_rbc = e2_right(end);
            e2_lbc = e2_left(1);
            e3_rbc = e3_right(end);
            e3_lbc = e3_left(1);
            
        end

        % Re-order and extend arrays to give reconstructions at all nx+1
        % interfaces, this uses periodicity since the nx+1st interface is the
        % same as the 1st interface.
        % "neg" means from below an interface, i.e., the reconstruction in the 
        % RHS of the cell to the interface's left.
        % "pos" means from above an interface, i.e., the reconstruction in the 
        % LHS of the cell to the interface's right.
        e1_pos = [e1_left; e1_rbc  ];
        e1_neg = [e1_lbc;  e1_right];
        e2_pos = [e2_left; e2_rbc  ];
        e2_neg = [e2_lbc;  e2_right];
        e3_pos = [e3_left; e3_rbc  ];
        e3_neg = [e3_lbc;  e3_right];
        
        
        rho_pos = obj.linearization_data.rho_pos{obj.linearization_data.t0_idx};
        rho_neg = obj.linearization_data.rho_neg{obj.linearization_data.t0_idx};
        u_pos   = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg   = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};
        E_pos   = obj.linearization_data.E_pos{obj.linearization_data.t0_idx};
        E_neg   = obj.linearization_data.E_neg{obj.linearization_data.t0_idx};

        % Compute physical flux at all interfaces
        [f1_pos, f2_pos, f3_pos, H_pos] = obj.flux_linearized(e1_pos, e2_pos, e3_pos, rho_pos, u_pos, E_pos);
        [f1_neg, f2_neg, f3_neg, H_neg] = obj.flux_linearized(e1_neg, e2_neg, e3_neg, rho_neg, u_neg, E_neg);

        % Compute Roe averages all all interfaces. See LeVeque p. 322--323.
        u_hat = (sqrt(rho_neg).*u_neg + sqrt(rho_pos).*u_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));
        H_hat = (sqrt(rho_neg).*H_neg + sqrt(rho_pos).*H_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));

        c_hat = sqrt( (obj.gamma-1)*( H_hat - 0.5*u_hat.^2 ) );

        % The jump for each conserved variable: 
        delta1 = e1_pos - e1_neg;
        delta2 = e2_pos - e2_neg;
        delta3 = e3_pos - e3_neg;

        % Coefficients representing jump in eigenvector basis.
        alpha2 = (obj.gamma - 1) .* ( (H_hat - u_hat.^2) .* delta1 + u_hat.*delta2 - delta3 ) ./ (c_hat.^2);

        % Based on LeVeque
        alpha3 = (delta2 + (c_hat - u_hat).*delta1 - c_hat.*alpha2) ./ (2*c_hat);
        alpha1 = delta1 - alpha2 - alpha3;

        % Smoothed absolute value, as in Harten's entropy fix (see LeVeque p. 326).
        lambda1_abs = obj.smoothed_absolute_value(u_hat - c_hat);
        lambda2_abs = obj.smoothed_absolute_value(u_hat);
        lambda3_abs = obj.smoothed_absolute_value(u_hat + c_hat);

        % Numerical flux.
        % See LeVeque, p.300--301 for eigenvalues and vectors.
        % Note the eigenvector components for q1 are just 1, so no need to
        % multiply by them.
        f_hat1 = 0.5*( f1_pos + f1_neg - ...
            (lambda1_abs.*alpha1 + ...
             lambda2_abs.*alpha2 + ...
             lambda3_abs.*alpha3 ) );
        % Eigenvector components for q2 are u-c, u, and u+c
        f_hat2 = 0.5*( f2_pos + f2_neg - ...
            (lambda1_abs.*alpha1.*(u_hat - c_hat) + ...
             lambda2_abs.*alpha2.*(u_hat ) + ...
             lambda3_abs.*alpha3.*(u_hat + c_hat)) );
        % Eigenvector components for q3 are H-u*c, 1/2*u^2, and H + u*c
        f_hat3 = 0.5*( f3_pos + f3_neg - ...
            (lambda1_abs.*alpha1.*(H_hat - u_hat.*c_hat) + ...
             lambda2_abs.*alpha2.*(0.5*u_hat.^2 ) + ...
             lambda3_abs.*alpha3.*(H_hat + u_hat.*c_hat)) );

        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat1(2:nx+1) - f_hat1(1:nx); ...     % First component.
               f_hat2(2:nx+1) - f_hat2(1:nx); ...     % Second component.
               f_hat3(2:nx+1) - f_hat3(1:nx)] / obj.mesh_pa.h; % Third component.
    end
    % End of spatial_discretization_linearized
    
    
    % See LeVeque eq. (14.43), p. 300.
    % We multiply the Jacobian matrix onto the vector [e1; e2; e3].
    function [f1, f2, f3, H] = flux_linearized(obj, e1, e2, e3, rho, u, E)
        
        gamma = obj.gamma;
        p = obj.pressure(rho, u, E);
        H = (E + p) ./ rho;
        
        f1 =                                                                e2;
        f2 = 0.5*(gamma - 3).*u.^2          .* e1 + (3 - gamma)*u        .* e2 + (gamma-1) .* e3;
        f3 = (0.5*(gamma - 1).*u.^3 - u.*H) .* e1 + (H - (gamma-1)*u.^2) .* e2 + gamma*u   .* e3;
    end
    
        
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
            
        end

        w_pos = [w_left; w_rbc];
        w_neg = [w_lbc;  w_right];
        
        rho_pos = obj.linearization_data.rho_pos{obj.linearization_data.t0_idx};
        rho_neg = obj.linearization_data.rho_neg{obj.linearization_data.t0_idx};
        u_pos   = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg   = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};
        E_pos   = obj.linearization_data.E_pos{obj.linearization_data.t0_idx};
        E_neg   = obj.linearization_data.E_neg{obj.linearization_data.t0_idx};

        % if t0_idx < obj.linearization_data.nt-1
        %     rho_pos = (2*obj.linearization_data.rho_pos{obj.linearization_data.t0_idx} + obj.linearization_data.rho_pos{obj.linearization_data.t0_idx+1})/3;
        %     rho_neg = (2*obj.linearization_data.rho_neg{obj.linearization_data.t0_idx} + obj.linearization_data.rho_neg{obj.linearization_data.t0_idx+1})/3;
        %     u_pos   = (2*obj.linearization_data.u_pos{obj.linearization_data.t0_idx}   + obj.linearization_data.u_pos{obj.linearization_data.t0_idx+1})/3;
        %     u_neg   = (2*obj.linearization_data.u_neg{obj.linearization_data.t0_idx}   + obj.linearization_data.u_neg{obj.linearization_data.t0_idx+1})/3;
        %     E_pos   = (2*obj.linearization_data.E_pos{obj.linearization_data.t0_idx}   + obj.linearization_data.E_pos{obj.linearization_data.t0_idx+1})/3;
        %     E_neg   = (2*obj.linearization_data.E_neg{obj.linearization_data.t0_idx}   + obj.linearization_data.E_neg{obj.linearization_data.t0_idx+1})/3;
        % end

        H_pos = (E_pos + obj.pressure(rho_pos, u_pos, E_pos)) ./ rho_pos;
        H_neg = (E_neg + obj.pressure(rho_neg, u_neg, E_neg)) ./ rho_neg;

        gamma = obj.gamma;
        
%         c_pos = sqrt( (gamma-1)*( H_pos - 0.5*u_pos.^2 ) );
%         c_neg = sqrt( (gamma-1)*( H_neg - 0.5*u_neg.^2 ) );
        
        c_pos = sqrt(gamma * obj.pressure(rho_pos, u_pos, E_pos) ./ rho_pos);
        c_neg = sqrt(gamma * obj.pressure(rho_neg, u_neg, E_neg) ./ rho_neg);
        
        assert(isreal(c_pos), 'c is not real...')
        
        % lambda1 = u - c
        if k == 1
            lambda_pos = u_pos - c_pos;
            lambda_neg = u_neg - c_neg;
        % lambda2 = u 
        elseif k == 2
            lambda_pos = u_pos;
            lambda_neg = u_neg;
        % lambda3 = u + c
        elseif k == 3
            lambda_pos = u_pos + c_pos;
            lambda_neg = u_neg + c_neg;
        else
            error('wave-speed index must be either 1 or 2 or 3.')
        end
        
        %% Compute the numerical flux
        % Compute physical flux at all interfaces
        f_pos = lambda_pos .* w_pos;
        f_neg = lambda_neg .* w_neg;

        % The jump for each conserved variable: 
        delta = w_pos - w_neg;

        %% Roe flux
        u_hat = (sqrt(rho_neg).*u_neg + sqrt(rho_pos).*u_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));
        H_hat = (sqrt(rho_neg).*H_neg + sqrt(rho_pos).*H_pos) ./ (sqrt(rho_neg) + sqrt(rho_pos));

        c_hat = sqrt( (obj.gamma-1)*( H_hat - 0.5*u_hat.^2 ) );

        % Smoothed absolute value, as in Harten's entropy fix (see LeVeque p. 326).
        if k == 1
            lambda_abs = obj.smoothed_absolute_value(u_hat - c_hat);
        elseif k == 2
            lambda_abs = obj.smoothed_absolute_value(u_hat);
        elseif k == 3
            lambda_abs = obj.smoothed_absolute_value(u_hat + c_hat);
        end
            
        f_hat = 0.5*( f_pos + f_neg - lambda_abs .* delta );

        % f_hat is a nx+1 dimensional vector with the numerical flux on all
        % nx+1 interfaces.
        % L(i) = - [f_hat(i+1/2) - f_hat(i-1/2)] / h.
        L = - [f_hat(2:nx+1) - f_hat(1:nx)] / obj.mesh_pa.h; 
    end
    % End of spatial_discretization_linearized_characteristic_variable
    
    % ic_description returns the name of the IC along with any parameters
    function [rho, u, E, ic_description] = initial_condition(obj, x)

        rho = zeros(size(x));
        u   = zeros(size(x));
        E   = zeros(size(x));
        gamma = obj.gamma;


        % Sod's problem.
        if strcmp(obj.pde_pa.ic_id, 'sod')
%             I = find(x < 0.5);
%             rho(I) = 1;
%             I = find(x >= 0.5);
%             rho(I) = 0.125;
% 
%             I = find(x < 0.5);
%             E(I) = 1/(gamma - 1);
%             I = find(x >= 0.5);
%             E(I) = 1/(gamma - 1) * 0.1;

            I = find(x < 0.5);
            rho(I) = 1;
            I = find(x >= 0.5);
            rho(I) = 1 - obj.pde_pa.ic_epsilon;

            I = find(x < 0.5);
            E(I) = 1/(gamma - 1);
            I = find(x >= 0.5);
            E(I) = 1/(gamma - 1) * (1 - obj.pde_pa.ic_epsilon);


            ic_description = sprintf('sod-eps%.2f', obj.pde_pa.ic_epsilon);
            

        % Shock-entropy problem. 
        elseif strcmp(obj.pde_pa.ic_id, 'shock-entropy')
            p = zeros(size(x));

            I = find(x < -4.0);
            rho(I) = 3.857143;
            u(I)   = 2.629369;
            p(I)   = 10.33333;

            I = find(x >= -4.0);
            rho(I) = 1 + 0.2*sin(pi*x(I));
            u(I)   = 0;
            p(I)   = 1;

            E = p/(gamma - 1) + 0.5*rho.*u.^2;
            
        elseif strcmp(obj.pde_pa.ic_id, 'idp1')
            beta = 5;
            rho = 1 + obj.pde_pa.ic_epsilon*exp(-beta * (x-2.5).^2);
            
            
            % change back to p = rho... and u = 0.

            % % if p is non-zero constant then get simple wave in density.
            % % in the simple wave, p and u are constant.
            % p   = 1*ones(size(x)) + 0.1*rho;
            % u = 0.2*p;

            p = rho;
            u = zeros(size(x));
            E = 1/(gamma - 1) * p + 0.5*rho .* u.^2;


            ic_description = sprintf('idp1-eps%.2f', obj.pde_pa.ic_epsilon);
        
        elseif strcmp(obj.pde_pa.ic_id, 'idp2')
            
            rho = 1 + obj.pde_pa.ic_epsilon*cos(pi*x/5);
            %rho = 0.5 + 0.001*cos(2*pi*x);
            
            p   = rho;
            
            u   = zeros(size(x));
            
            E = 1/(gamma - 1) * p + 0.5*rho .* u.^2;
            
        end

    end
    % End of initial condition
    
    % See LeVeque eq. (14.39), p. 299.
    function lambda = wave_speeds(obj, q)
        nx = obj.mesh_pa.nx;
        gamma = obj.gamma;
        
        rho = q(1:nx);
        u   = q(nx+1:2*nx) ./ rho;
        E   = q(2*nx+1:3*nx);
        
        c = sqrt( gamma .* obj.pressure(rho, u, E) ./ rho );
        
        lambda = [u - c; u; u + c];
    end
    
    
    % Map from characteristic variables b to primitive variables a. This
    % is achieved by applying R to b.
    % q is the state used to evaluate the eigenvectors.
    % See LeVeque eq. (14.46), p. 301.
    function a = right_eigenvector_map(obj, q, b) 
        nx = obj.mesh_pa.nx;
        gamma = obj.gamma;
        
        rho = q(1:nx);
        u   = q(nx+1:2*nx) ./ rho;
        E   = q(2*nx+1:3*nx);
        
        p = obj.pressure(rho, u, E);
        H = (E + p) ./ rho;
        c = sqrt( gamma .* p ./ rho );
        
        b1 = b(1:nx);
        b2 = b(nx+1:2*nx);
        b3 = b(2*nx+1:3*nx);
        
        % First row of R multiplied with b
        a1 =               b1 +             b2 +               b3; 
        % Second row of R multiplied with b
        a2 = (u-c)      .* b1 + u        .* b2 + (u+c)      .* b3;
        % Third row of R multiplied with b
        a3 = (H - u.*c) .* b1 + 0.5*u.^2 .* b2 + (H + u.*c) .* b3;
        
        a = [a1; a2; a3];
    end
    
    % Map from primitive variables a to characteristic variables b. This
    % is achieved by applying R^-1 to a.
    % q is the state used to evaluate the eigenvectors.
    function b = left_eigenvector_map(obj, q, a) 
        nx = obj.mesh_pa.nx;
        gamma = obj.gamma;
        
        rho = q(1:nx);
        u   = q(nx+1:2*nx) ./ rho;
        E   = q(2*nx+1:3*nx);
        
        p = obj.pressure(rho, u, E);
        H = (E + p) ./ rho;
        c = sqrt( gamma .* p ./ rho );
        
        a1 = a(1:nx);
        a2 = a(nx+1:2*nx);
        a3 = a(2*nx+1:3*nx);
        
        g = 1/(gamma - 1);
        
        % First row of R multiplied with b
        b1 = (H + g*c.*(u - c) ) .* a1 + (-u - g*c) .* a2 +   a3; 
        b1 = b1 * 0.5*(gamma - 1) ./ (c.^2);
        % Second row of R multiplied with b
        b2 = (4*g*c.^2 - 2*H)    .* a1 + 2*u        .* a2 - 2*a3;
        b2 = b2 * 0.5*(gamma - 1) ./ (c.^2);
        % Third row of R multiplied with b
        b3 = (H - g*c.*(u + c) ) .* a1 + (-u + g*c) .* a2 +   a3; % Typo in formula in (3,1) in place I got this from...
        b3 = b3 * 0.5*(gamma - 1) ./ (c.^2);
        
        b = [b1; b2; b3];     
    end

    % Return the kth linearized wave-speed based at time t0 based on the
    % stored linearization data.
    % "neg" means from below an interface, i.e., the reconstruction in the 
    % RHS of the cell to the interface's left.
    % "pos" means from above an interface, i.e., the reconstruction in the 
    % LHS of the cell to the interface's right.
    function [lambda_k_neg, lambda_k_pos] = linearized_wave_speed(obj, k) 
        assert(~isempty(obj.linearization_data), 'Cannot compute linearized wave-speed without linearization data');

        rho_pos = obj.linearization_data.rho_pos{obj.linearization_data.t0_idx};
        rho_neg = obj.linearization_data.rho_neg{obj.linearization_data.t0_idx};
        u_pos   = obj.linearization_data.u_pos{obj.linearization_data.t0_idx};
        u_neg   = obj.linearization_data.u_neg{obj.linearization_data.t0_idx};
        E_pos   = obj.linearization_data.E_pos{obj.linearization_data.t0_idx};
        E_neg   = obj.linearization_data.E_neg{obj.linearization_data.t0_idx};
        c_pos = sqrt(obj.gamma * obj.pressure(rho_pos, u_pos, E_pos) ./ rho_pos);
        c_neg = sqrt(obj.gamma * obj.pressure(rho_neg, u_neg, E_neg) ./ rho_neg);

        if k == 1
            lambda_k_neg = u_neg - c_neg;
            lambda_k_pos = u_pos - c_pos;
        elseif k == 2
            lambda_k_neg = u_neg;
            lambda_k_pos = u_pos;
        elseif k == 3
            lambda_k_neg = u_neg + c_neg;
            lambda_k_pos = u_pos + c_pos;
        end
        
    end
    
    end
    % End of methods

end


% %         b = randn(my_cons_law.m*nx, 1);
% %         Rb = my_cons_law.right_eigenvector_map(q(:, 1), b);
% %         a = my_cons_law.left_eigenvector_map( q(:, 1), Rb);
% %         norm(b - a)
%         rho = rand(1);
%         u   = rand(1);
%         E   = rand(1);
% 
%         g = 1/(gamma - 1);
%         p = obj.pressure(rho, u, E);
%         H = (E + p) ./ rho;
%         c = sqrt( gamma .* p ./ rho );
% 
%         a1 = 1;
%         a2 = 1;
%         a3 = 1;
%         Rinv = 0.5*(gamma-1)/(c.^2).*[...
%             (H + g*c.*(u - c) ) .* a1, + (-u - g*c) .* a2, +   a3; ...
%             (4*g*c.^2 - 2*H)    .* a1, + 2*u        .* a2, - 2*a3; ...
%             (H - g*c.*(u + c) ) .* a1, + (-u + g*c) .* a2, +   a3;
%         ];
% 
% 
%         b1 = 1;
%         b2 = 1;
%         b3 = 1;
%         R = [...
%                             b1, +             b2, +               b3; ...
%             (u-c)      .* b1, + u        .* b2, + (u+c)      .* b3; ...
%          (H - u.*c) .* b1, + 0.5*u.^2 .* b2, + (H + u.*c) .* b3];
% 
%         norm(Rinv*R - eye(3), inf)
% 

