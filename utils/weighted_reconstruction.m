% Class defining weighted reconstructions of 1D functions.
%
% TODO: Set up some kind of way for the user to specify which WENO weights
% they want and which parameters... Also set a flag for the gradients to be
% returned only if the classical weights are used...
classdef weighted_reconstruction < handle
    
    properties
        k;
        p;
        nx;
        id;

        d_left;
        d_right;
        B;
        reconstruct_matrices_left;
        reconstruct_matrices_right;
    end
    
    methods
        
        % The methods implemented below are:
        % weighted_reconstruction
        % WENO_weights
        % smoothness_indicators
        % interface_reconstruction
        % shifted_reconstruction_matrices
        
        % Reconstruct solution at left and right interace in every cell. 
        % This can be done in a standard way using either linear weights or
        % WENO weights. Otherwise an alternative set of weights can be
        % parsed and these used instead.
        % The WENO weights can also be returned.
        function [u_left_interface, u_right_interface, weights] ...
                    = interface_reconstruction(obj, u_bar, ...
                                    weights_left, weights_right)
                
            weights = struct();
                                
            % Do reconstruction with the non-standard weights if parsed
            if nargin == 2
                standard_reconstruction = true;
            elseif nargin == 4
                standard_reconstruction = false;
            end
                                
            % Inside each cell we reconstruct u at left- and right-hand interfaces,
            % and we do this on k different stencils. 
            u_left_interface  = zeros(obj.nx, obj.k);
            u_right_interface = zeros(obj.nx, obj.k);

            % Do the k different reconstructions, each with a different left shift
            % ell.
            for ell = 0:obj.k-1
                u_left_interface(:, ell+1)  = obj.reconstruct_matrices_left{ell+1}  * u_bar;
                u_right_interface(:, ell+1) = obj.reconstruct_matrices_right{ell+1} * u_bar;
            end

            % Combine the k different k cell reconstructions to create one big 2k-1 
            % cell reconstruction.
            if standard_reconstruction
                
                if strcmp(obj.id, 'linear')
                    % The sums here are computed by inner products
                    u_left_interface  = u_left_interface  * obj.d_left; 
                    u_right_interface = u_right_interface * obj.d_right;

                elseif strcmp(obj.id, 'WENO')
                    [omega_left, omega_right] = obj.WENO_weights(u_bar);
                    u_left_interface  = sum(u_left_interface  .* omega_left, 2);
                    u_right_interface = sum(u_right_interface .* omega_right, 2);

                    if nargout == 3
                        weights.left  = omega_left;
                        weights.right = omega_right;
                    end
                end
            
                % Reconstruction using the parsed weights.     
            elseif ~standard_reconstruction
                u_left_interface  = sum(u_left_interface  .* weights_left, 2);
                u_right_interface = sum(u_right_interface .* weights_right, 2);
                
            end

        end
        % End of: interface_reconstruction
        
        % Reconstruct error e at left and right interace in every cell using
        % reconstruction linearized about u_bar.
        function [e_left_interface, e_right_interface] ...
            = interface_reconstruction_WENO_linearized(obj, u_bar, e_bar)


            % Inside each cell we reconstruct u and e at left- and right-hand 
            % interfaces, and we do this on k different stencils. 
            u_left_interface  = zeros(obj.nx, obj.k);
            u_right_interface = zeros(obj.nx, obj.k);

            e_left_interface  = zeros(obj.nx, obj.k);
            e_right_interface = zeros(obj.nx, obj.k);

            % Do the k different reconstructions, each with a different left shift
            % ell. Do this for both u and e.
            for ell = 0:obj.k-1
                u_left_interface(:, ell+1)  = obj.reconstruct_matrices_left{ell+1}  * u_bar;
                u_right_interface(:, ell+1) = obj.reconstruct_matrices_right{ell+1} * u_bar;

                e_left_interface(:, ell+1)  = obj.reconstruct_matrices_left{ell+1}  * e_bar;
                e_right_interface(:, ell+1) = obj.reconstruct_matrices_right{ell+1} * e_bar;
            end

            % Compute WENO weights based on u, and their gradients in the direction
            % of e_bar.
            [omega_left,      omega_right, ...
             omega_left_grad, omega_right_grad] ...
                = obj.WENO_weights(u_bar, e_bar);

            % Combine the k different k cell reconstructions to create one big 2k-1
            % cell reconstruction.

            % The linearized reconstruction has two components. The first fixes the
            % u_bar-based WENO weights and reconstructions e with them, and the
            % second does a reconstruction with u based on the gradient of the WENO
            % weights in the direcion of e_bar
            e_left_interface  = sum(e_left_interface  .* omega_left, 2)  + sum(u_left_interface  .* omega_left_grad, 2);
            e_right_interface = sum(e_right_interface .* omega_right, 2) + sum(u_right_interface  .* omega_right_grad, 2);

        end
        % End of: interface_reconstruction_WENO_linearized
        
        % Constructor
        function obj = weighted_reconstruction(p_, nx_, id_)
            
            
            assert(strcmp(id_, 'linear') || ...
                   strcmp(id_, 'WENO'), ...
                   '''id'' must be linear or WENO');
            
            assert(~(strcmp(id_, 'WENO') && p_ == 1), ...
                '1st-order WENO reconstruction does not exist');
               
            obj.id = id_;
            
            obj.p  = p_;
            obj.k = (p_ + 1)/2;
            obj.nx = nx_;

            assert(p_ > 0, 'require p > 0')
            assert(p_ <= 5, 'p > 5 not implemented')

            % Get linear weights
            [obj.d_left, obj.d_right] = obj.linear_weights();
            
            % Get smoothness indicator matrices
            obj.B = obj.smoothness_stencil_matrices();
            
            % Get matrices for interface reconstructions.
            [obj.reconstruct_matrices_left, ...
             obj.reconstruct_matrices_right] = ...
                                    obj.shifted_reconstruction_matrices();
            
        end
        % End of: weighted_reconstruction
         
        
        % See Shu (1997) eqns. (2.58) & (2.59), page 18, and also page. 20 
        % for the left interface weights. 
        %
        % If e_bar is given, then also the gradient of the WENO weights at 
        % u_bar in the direction of e_bar is given.
        function [omega_left,      omega_right, ...
                  omega_left_grad, omega_right_grad] = ...
                    WENO_weights(obj, u_bar, e_bar)

            % Do we need to compute gradients also?
            if nargin > 2; compute_gradients = true; else; compute_gradients = false; end
               
            % Get smoothness measures
            if ~compute_gradients
                beta = obj.smoothness_indicators(u_bar);
            else
                [beta, beta_gradient] = obj.smoothness_indicators(u_bar, e_bar);
            end
            

            epsilon = 1e-6;
            %epsilon = 1e-12;
            %epsilon = pa.h^1;

            pp = 1;

            % From Hesthaven (2017), p. 338: 
            % typically " = 10^-6. The appropriate value of " is closely related to 
            % the magnitude of beta and, for very accurate solutions, a smaller 
            % value may be needed
            % The value of p can be changed to modify the impact of the smoothness 
            % indicator, although a value of p = 1 is often used.
            
            alpha_left   = obj.d_left.'  ./ (epsilon + beta).^(2*pp); 
            alpha_right  = obj.d_right.' ./ (epsilon + beta).^(2*pp);        
        
            
%             epsilon = 1e-6;
%             pp = 1; % Scheme is then WENOZ (see Sec. 2.3.2, p. 3435)
%             pp = 3/4; % They seem to make some conclusion that maybe pp = 0.75 is
%             % best...
%             tau = abs(beta(:, 1) - beta(:,end));
%             alpha_left   = obj.d_left.'  .* (1 + (tau ./ (epsilon + beta).^pp) ); 
%             alpha_right  = obj.d_right.' .* (1 + (tau ./ (epsilon + beta).^pp) ); 

            
            omega_left   = alpha_left  ./ sum(alpha_left, 2);
            omega_right  = alpha_right ./ sum(alpha_right, 2);
        
            
            % Evaluate gradients 
            if compute_gradients
               
                % Form the sums that appear in all left and right gradients.
                sum_left  = sum( alpha_left  .* sqrt(alpha_left  ./ obj.d_left.')  .* beta_gradient, 2);
                sum_right = sum( alpha_right .* sqrt(alpha_right ./ obj.d_right.') .* beta_gradient, 2);

                omega_left_grad  = 2 * omega_left ./ alpha_left .* ( ...
                      - alpha_left  .* sqrt(alpha_left  ./ obj.d_left.')  .* beta_gradient ...
                      + omega_left  .* sum_left );

                omega_right_grad = 2 * omega_right ./ alpha_right .* ( ...
                      - alpha_right .* sqrt(alpha_right ./ obj.d_right.') .* beta_gradient ...
                      + omega_right .* sum_right );


            %         I = find(sum(beta, 2) < 1e-5);
            %         [size(beta, 1), numel(I), min(sum(beta, 2)), median(sum(beta, 2)), max(sum(beta, 2))]
            %     
            %         figure(1)
            %         clf
            %         plot(omega_left_grad(:, 1), '-bo')
            %         figure(2)
            %         clf
            %         %plot(sum(beta, 2), '-rx')
            %         %hold on
            %         plot(sum(beta_gradient, 2), '-c>')
            %         
            %         pause
            %         
            % %         plot(omega_left_grad, '-bo')
            % %         hold on
            %         omega_left_grad(I, :)  = 0.0;
            %         omega_right_grad(I, :) = 0.0;
            % %         plot(omega_left_grad, 'rx')
                
            end
            

            % Some exprimenting with non-classical WENO weights.
            % tau = abs(beta(:, 1) - beta(:,end));
            % pp = 3/4;
            % alpha_left   = obj.d_left.'  .* (1 + (tau ./ (epsilon + beta)).^pp ); 
            % alpha_right  = obj.d_right.' .* (1 + (tau ./ (epsilon + beta)).^pp ); 
            % omega_left   = alpha_left  ./ sum(alpha_left, 2);
            % omega_right  = alpha_right ./ sum(alpha_right, 2);
            % 
            % epsilon = 1e-6;
            % %epsilon = pa.h;
            % tau = abs(beta(:, 1) - beta(:,end));
            % alpha_left   = obj.d_left.'  .* (1 + (tau ./ (epsilon + beta)) ); 
            % alpha_right  = obj.d_right.' .* (1 + (tau ./ (epsilon + beta)) ); 
            % omega_left   = alpha_left  ./ sum(alpha_left, 2);
            % omega_right  = alpha_right ./ sum(alpha_right, 2);


            % Here is a non-vectorized way to compute the above weights.        
            %     alpha_left  = zeros(pa.nx, k);
            %     alpha_right = zeros(pa.nx, k);
            %     alpha_left_sum  = 0;
            %     alpha_right_sum = 0;
            %     for ell = 1:k
            %        alpha_left(:, ell)  = d_left(ell)  ./ (epsilon + beta(:, ell)).^2; 
            %        alpha_right(:, ell) = d_right(ell) ./ (epsilon + beta(:, ell)).^2; 
            %        
            %        alpha_left_sum  = alpha_left_sum + alpha_left(:, ell);
            %        alpha_right_sum = alpha_right_sum + alpha_right(:, ell);
            %     end
            %     omega_left  = alpha_left  ./ alpha_left_sum;
            %     omega_right = alpha_right ./ alpha_right_sum;

        end
        
        
        % Evaluate the smoothness indicators. (See Shu (1997) eqns (2.62) and
        % (2.63)).
        %
        % If e_bar is given, then the gradient of the smoothness measures
        % at u_bar in the direction of e_bar is computed also.
        function [beta, beta_gradient] = smoothness_indicators(obj, u_bar, e_bar)
            if obj.k == 1
                beta = zeros(size(u_bar));
                
            elseif obj.k == 2
                beta = [(obj.B{1}*u_bar).^2, (obj.B{2}*u_bar).^2];
                
            elseif obj.k == 3
                beta = [13/12*(obj.B{1,1}*u_bar).^2 + 1/4*(obj.B{1,2}*u_bar).^2, ...
                        13/12*(obj.B{2,1}*u_bar).^2 + 1/4*(obj.B{2,2}*u_bar).^2, ...
                        13/12*(obj.B{3,1}*u_bar).^2 + 1/4*(obj.B{3,2}*u_bar).^2];
            end
            
            % Compute gradients
            if nargin > 2
                if obj.k == 1
                    beta_gradient = zeros(size(u_bar));

                elseif obj.k == 2
                    beta_gradient = [2*(obj.B{1}*u_bar).*(obj.B{1}*e_bar), 2*(obj.B{2}*u_bar).*(obj.B{2}*e_bar)];

                elseif obj.k == 3
                   beta_gradient = [2*13/12*(obj.B{1,1}*u_bar).*(pa.B{1,1}*e_bar) + 2*1/4*(obj.B{1,2}*u_bar).*(obj.B{1,2}*e_bar), ...
                                    2*13/12*(obj.B{2,1}*u_bar).*(pa.B{2,1}*e_bar) + 2*1/4*(obj.B{2,2}*u_bar).*(obj.B{2,2}*e_bar), ...
                                    2*13/12*(obj.B{3,1}*u_bar).*(pa.B{3,1}*e_bar) + 2*1/4*(obj.B{3,2}*u_bar).*(obj.B{3,2}*e_bar)];     
                end
            end
            
        end
        % End of: smoothness_indicators
        
        
        

    
        % Linear weights for combining k different reconstructions.
        % See Shu (1997), p. 17. The left interface coefficients are just the
        % reverse of the left ones (see p. 20).
        function [d_left, d_right] = linear_weights(obj)
            if obj.k == 1
                d_right = 1;
            elseif obj.k == 2
                d_right = [2/3; 1/3];
            elseif obj.k == 3
                d_right = [3/10; 3/5; 1/10];
            end

            d_left = flipud(d_right);
        end
    
        
        % Build matrices that when applied to a vector of cell avaerages of u will 
        % reconstruct u at both left and right interfaces in every cell. There are 
        % k such matrices since there are k such reconstructions.
        function [reconstruct_matrices_left_interface, ...
                  reconstruct_matrices_right_interface] = ...
                                       shifted_reconstruction_matrices(obj)

            spatial_problem          = struct();
            spatial_problem.nDOFs    = obj.nx;
            spatial_problem.BCs.type = 'periodic';

            reconstruct_matrices_left_interface  = cell(obj.k, 1);
            reconstruct_matrices_right_interface = cell(obj.k, 1);

            % Left-shift of stencil. Start from the right-most stencil which uses
            % ell = 0. There are k different stencils
            for ell = 0:obj.k-1

                % Get left- and right-interface reconstruction weights.
                [c_left, c_right] = obj.interface_reconstruction_weights(ell);

                % Corresponding right-shift of stencil given that the stencil contains k cells.
                r = (obj.k - 1) - ell; 

                % Nodes used in the stencil.
                reconstruct_nodes = (-ell:r)';

                % Left interface
                spatial_problem.stencil                     = {{c_left;  reconstruct_nodes}};             
                reconstruct_matrices_left_interface{ell+1}  = get_toeplitz_spatial_disc(spatial_problem);

                % Right interface
                spatial_problem.stencil                     = {{c_right; reconstruct_nodes}};             
                reconstruct_matrices_right_interface{ell+1} = get_toeplitz_spatial_disc(spatial_problem); 

            end  
        end
        % End of: shifted_reconstruction_matrices
        
        
        % Get weights for a reconstruction at left and right cell interfaces that 
        % uses a stencil with k cells and is shifted ell cells to the left. 
        % These weights are just Lagrange polynomials evaluated at cell interfaces.
        % See Shu (1997), Table 2.1, page 8.
        function [c_left, c_right] = interface_reconstruction_weights(obj, ell)

            if ell < 0 || ell > obj.k-1
                error('Left shift ell of stencil needs to be in [0,k-1]')
            end

            c_right = obj.right_interface_reconstruction_weights(ell);

            % The stencil coefficients for the left-interface reconstruction are 
            % the same as for the right interface, but with the value of ell one
            % less. And for the ell = 0 case, the coefficients are given by the
            % reverse of the right interface values from the ell=k-1 case.
            if ell > 0
                c_left = obj.right_interface_reconstruction_weights(ell-1);

            % The output of right_interface_reconstruction_weights() is a column
            % vector, so do flipud to reverse the ordering.
            elseif ell == 0

                c_left = flipud( obj.right_interface_reconstruction_weights(obj.k-1) );
            end

        end
        % End of: interface_reconstruction_weights
        
        
        % The left-shift index ell goes from 0 to k-1. I ignore here the ell=-1
        % case because this is just shown in Shu's table because it relates to
        % how the left-interface reconstruction coefficients are computed.
        function c_right = right_interface_reconstruction_weights(obj, ell)

            assert(ell >= 0 && ell <= obj.k-1, 'Left shift ell needs to be in [0,k-1]')

            if obj.k == 1
                c_right = 1;

            elseif obj.k == 2
                if ell == 0
                    c_right = [ 1/2;  1/2];

                elseif ell == 1
                    c_right = [-1/2;  3/2];     
                end

            elseif obj.k == 3
                if ell == 0
                    c_right = [ 1/3;  5/6; -1/6];

                elseif ell == 1
                    c_right = [-1/6;  5/6;  1/3]; 

                elseif ell == 2
                    c_right = [ 1/3; -7/6;  11/6];   
                end
                
           elseif obj.k == 4
                if ell == 0
                    c_right = [ 1/4;  13/12;  -5/12;  1/12];

                elseif ell == 1
                    c_right = [-1/12;  7/12;   7/12; -1/12]; 

                elseif ell == 2
                    c_right = [ 1/12; -5/12;  13/12;  1/4];
                    
                elseif ell == 3
                    c_right = [-1/4;   13/12; -23/12; 25/12];
                end
                
            elseif obj.k == 5
                if ell == 0
                    c_right = [  1/5;  77/60; -43/60;  17/60; -1/20];

                elseif ell == 1
                    c_right = [-1/20;  9/20;   47/60; -13/60;  1/30]; 

                elseif ell == 2
                    c_right = [ 1/30; -13/60;  47/60;   9/20; -1/20];
                    
                elseif ell == 3
                    c_right = [-1/20; 17/60;  -43/60;  77/60;  1/5];
                    
                elseif ell == 4
                    c_right = [ 1/5; -21/20;  137/60; -163/60; 137/60];
                end    

            end
        end
        % End of: right_interface_reconstruction_weights


        % Get matrices than implement the stencils used to evaluate smoothness
        % measures (in this way, we don't have to worry about splicing the solution 
        % vector at boundaries, etc. since it's all taken care of by a MATVEC).
        % After performing the MATVEC we still need to square the resulting vector
        % and weight it with the appropriate coefficients (this is all done in the
        % smoothness_indicators function).
        function B = smoothness_stencil_matrices(obj)

            spatial_problem.nDOFs    = obj.nx;
            spatial_problem.BCs.type = 'periodic';

            % There is no smoothness measure for 1st-order reconstructions.
            if obj.k == 1
                B = [];

            elseif obj.k == 2
                B       = cell(2, 1);

                %% First smoothness measure, ell=0
                nodes   = [ 0; 1];
                weights = [-1; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{1}    = get_toeplitz_spatial_disc(spatial_problem); 

                %% Second smoothness measure, ell=1
                nodes   = [-1; 0];
                weights = [-1; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{2}    = get_toeplitz_spatial_disc(spatial_problem); 

            elseif obj.k == 3
                B = cell(3, 2);

                %% First smoothness measure, ell=0
                % First term
                nodes   = [0;  1; 2];
                weights = [1; -2; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{1,1}  = get_toeplitz_spatial_disc(spatial_problem); 

                % Second term
                nodes   = [0;  1; 2];
                weights = [3; -4; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{1,2}  = get_toeplitz_spatial_disc(spatial_problem); 

                %% Second smoothness measure, ell=1
                % First term
                nodes   = [-1;  0; 1];
                weights = [ 1; -2; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{2,1}  = get_toeplitz_spatial_disc(spatial_problem); 

                % Second term
                nodes   = [-1; 0; 1];
                weights = [ 1; 0; -1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{2,2}  = get_toeplitz_spatial_disc(spatial_problem);

                %% Third smoothness measure, ell=2
                % First term
                nodes   = [-2; -1; 0];
                weights = [ 1; -2; 1];
                spatial_problem.stencil = {{weights; nodes}};             
                B{3,1}  = get_toeplitz_spatial_disc(spatial_problem); 

                % Second term
                nodes   = [-2; -1; 0];
                weights = [ 1; -4; 3];
                spatial_problem.stencil = {{weights; nodes}};             
                B{3,2}  = get_toeplitz_spatial_disc(spatial_problem);

            end
        end
        % End of: smoothness_stencil_matrices

    end
    % End of: methods
    
end
% End of: file