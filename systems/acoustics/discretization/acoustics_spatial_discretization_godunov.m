% Implements the spatial discretization for the acoustic system that arises
% from solving exactly the Riemann problem at each cell interface.
%
% The wave-speed and impedance at cell centers are expected in the struct 
% disc_pa stored as fields c_cc and Z_cc, respectively. 
%
% q should be a stacked vector with the first nx components the cell
% averages of the pressure and the second nx the cell averages of the
% velocity. 
%
% Uses periodic boundaries, but is not very hard to adapt to implement
% inflow/outflow boundaries. 
function L = acoustics_spatial_discretization_godunov(q, mesh_pa, disc_pa)
    
    nx = mesh_pa.nx;
    
    % Extract q1 and q2 from stacked q vector.
    p = q(1:nx);
    u = q(nx+1:2*nx);

    % Reconstructions at left and right interfaces of all nx cells. 
    % Reconstructions are 1st-order: The solution is assumed constant in 
    % each cell. 
    p_left  = p;
    p_right = p;
    u_left  = u;
    u_right = u;
    
    % Re-order and extend arrays to give reconstructions on both sides of 
    % all nx+1 interfaces.
    % This uses periodicity since the nx+1st interface is the same as the 
    % 1st interface. "neg" means from below an interface, i.e., the 
    % reconstruction in the RHS of the cell to the interface's left.
    % "pos" means from above an interface, i.e., the reconstruction in the 
    % LHS of the cell to the interface's right.
    p_pos = [p_left;      p_left(1)]; % p on all nx+1 interfaces from ABOVE
    p_neg = [p_right(nx); p_right  ]; % p on all nx+1 interfaces from BELOW
    u_pos = [u_left;      u_left(1)]; % u on all nx+1 interfaces from ABOVE
    u_neg = [u_right(nx); u_right  ]; % u on all nx+1 interfaces from BELOW


    % Similarly for the material parameters, store them above and below
    % every interface. 
    % Impedance
    Z_left  = disc_pa.Z_cc;
    Z_right = disc_pa.Z_cc;
    Z_pos   = [Z_left;      Z_left(1)];
    Z_neg   = [Z_right(nx); Z_right  ];
    % Sound speed
    c_left  = disc_pa.c_cc;
    c_right = disc_pa.c_cc;
    c_pos   = [c_left;      c_left(1)];
    c_neg   = [c_right(nx); c_right  ];

    % Strength of waves.
    alpha1 = ( -(p_pos - p_neg) + Z_pos .* (u_pos - u_neg) ) ./ (Z_neg + Z_pos);
    alpha2 = (  (p_pos - p_neg) + Z_neg .* (u_pos - u_neg) ) ./ (Z_neg + Z_pos);
    
    %[norm(alpha1), norm(alpha2)]

    W1 = [alpha1 .* -Z_neg, alpha1]; % LEFT-going waves emanating from all nx+1 interfaces
    W2 = [alpha2 .*  Z_pos, alpha2]; % RIGHT-going waves emanating from all nx+1 interfaces
    % (note that we don't actually need the LEFT-going wave eminating from
    % the first interface nor the right-going wave emanating from the
    % nx+1st interface, but for simplicity of the implementation we compute
    % these anyway).
    
    % The flux in each cell if the sum of the right-going wave entering it
    % from the LHS interface of that cell, and the left-going wave entering
    % it from the RHS interface of that cell. Each is multiplied by the
    % speed of the wave. 
    L = -[c_pos(1:nx) .* W2(1:nx, 1) - c_neg(2:nx+1) .* W1(2:nx+1, 1); ...          % First component.
          c_pos(1:nx) .* W2(1:nx, 2) - c_neg(2:nx+1) .* W1(2:nx+1, 2) ] / mesh_pa.h; % Second component.
      
    % Add high-resolution corrections to the Godunov discretization.
    % Implements formula (9.67) for all nx+1 interfaces. Noting that |s^1| =
    % c_neg and |s^2| = c_pos (c is always positive).
    if isfield(disc_pa, 'high_res') && disc_pa.high_res
        h  = mesh_pa.h;
        dt = mesh_pa.dt;
        
        % First component of F_tilde
        F_tilde1 = 0.5*c_pos(1:nx+1) .* (1 - dt/h * c_pos(1:nx+1)) .* W2(1:nx+1, 1) + 0.5*c_neg(1:nx+1) .* (1 - dt/h * c_neg(1:nx+1)) .* W1(1:nx+1, 1);
        % Second component of F_tilde
        F_tilde2 = 0.5*c_pos(1:nx+1) .* (1 - dt/h * c_pos(1:nx+1)) .* W2(1:nx+1, 2) + 0.5*c_neg(1:nx+1) .* (1 - dt/h * c_neg(1:nx+1)) .* W1(1:nx+1, 2);

        % Interestingly, this correction is conservative (while the main term
        % is not)... This is probably related to the Taylor series discussion
        % following (9.68).
        L = L - [F_tilde1(2:nx+1) - F_tilde1(1:nx); F_tilde2(2:nx+1) - F_tilde2(1:nx)] / h;
    end

end