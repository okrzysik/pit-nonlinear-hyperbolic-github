% Build time-stepping operators for Godunov discretizations of linear
% advection equations wave wave-speeds \pm c.
%
% Alternatively, use LW (Lax--Wendroff) discretizations of 
% linear advection equations wave wave-speeds \pm c: See eq. (6.4) from 
% LeVeque (2004).

function [Phi_neg, Phi_pos] = Phi_blocks_advection_godunov(mesh_pa, disc_pa)

    % Discretize with LW or Godunov?
    LW = ~true;
    if isfield(disc_pa, 'high_res') && disc_pa.high_res
        LW = true;
    end

    nx = mesh_pa.nx;
    h  = mesh_pa.h;
    dt = mesh_pa.dt;

    % Build and the spatial discretization operators.
    I  = speye(nx);

    % Godunov disc.
    if ~LW
        % Build D1: First-order upwind derivative
        derivative = 1;
        bias = 'U';
        approx_order = 1;
        [p_nodes, p_weights] = FD_stencil(derivative, bias, approx_order);
        p_weights = p_weights/(h^derivative);
        spatial_problem.stencil  = {{p_weights; p_nodes}};             
        spatial_problem.nDOFs = nx;
        spatial_problem.BCs.type = 'periodic';
        L = get_toeplitz_spatial_disc(spatial_problem); 

        % Assemble discretizations as matrices.
        Phi_neg = I - dt * (-disc_pa.c_cc) .* (-L.');
        Phi_pos = I - dt * (+disc_pa.c_cc) .* (+L);
        
    % LW disc.
    else
    
        % Build and the spatial discretization operators.
        I  = speye(nx);

        % Build L: Central second-order first derivative
        L_nodes   = [-1; 0; 1];
        L_weights = [-1; 0; 1];
        spatial_problem.stencil  = {{L_weights; L_nodes}};             
        spatial_problem.nDOFs = nx;
        spatial_problem.BCs.type = 'periodic';
        L = get_toeplitz_spatial_disc(spatial_problem); 

        % Build D: Central second-order second derivative
        D_nodes   = [-1;  0; 1];
        D_weights = [ 1; -2; 1];
        spatial_problem.stencil  = {{D_weights; D_nodes}};             
        spatial_problem.nDOFs = nx;
        spatial_problem.BCs.type = 'periodic';
        D = get_toeplitz_spatial_disc(spatial_problem); 

        % Assemble discretizations as matrices.

        Phi_neg = I - 0.5 * (dt / h * -disc_pa.c_cc) .* L + 0.5 * (dt / h * -disc_pa.c_cc).^2 .* D;
        Phi_pos = I - 0.5 * (dt / h * +disc_pa.c_cc) .* L + 0.5 * (dt / h * +disc_pa.c_cc).^2 .* D;
        
    end

end
