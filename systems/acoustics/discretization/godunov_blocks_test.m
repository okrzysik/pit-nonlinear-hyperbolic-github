% This is a check to make sure I did the algebra right in transforming the
% discretization from primitive varibles into characteristic variables.
% Everything seemig OK...

clc
clear


CFL = 0.85;

nx = 2.^8;
h = 1/nx;

dt = CFL * h;

x = linspace(0, 1, nx);


pde_pa.mat_param_id = 0;
pde_pa.mat_param_id = 1;
pde_pa.mat_param_id = 2;
%pde_pa.mat_param_id = 3;
%pde_pa.mat_param_id = 4;
%pde_pa.mat_param_id = 5;
pde_pa.mat_param_id = 6;

pde_pa.mat_param_num_layers = 16;

[c, Z] = acoustics_material_parameters(x, pde_pa);

[Phi_pp, Phi_pu, Phi_up, Phi_uu] = Phi_blocks_acoustics_godunov(dt, h, c, Z, nx);

[Phi_hat_11, Phi_hat_12, Phi_hat_21, Phi_hat_22] = Phi_hat_blocks_acoustics_godunov(dt, h, c, Z, nx);


Phi      = [Phi_pp, Phi_pu; Phi_up, Phi_uu];
Phi_hat = [Phi_hat_11, Phi_hat_12; Phi_hat_21, Phi_hat_22];

Z_mat     = diag(Z);
Z_mat_inv = diag(1./Z);
I    = speye(nx);
R    =     [-Z_mat, Z_mat; I, I];
Rinv = 0.5*[-Z_mat_inv, I; Z_mat_inv, I];

Delta = Rinv * Phi * R - Phi_hat;

norm(Delta, inf)

%spy(Delta)

Phi_11_hat_true = 0.5*( (Phi_uu + Z_mat_inv * Phi_pp * Z_mat) - (Phi_up*Z_mat + Z_mat_inv*Phi_pu) );
Phi_12_hat_true = 0.5*( (Phi_uu - Z_mat_inv * Phi_pp * Z_mat) + (Phi_up*Z_mat - Z_mat_inv*Phi_pu) );
Phi_21_hat_true = 0.5*( (Phi_uu - Z_mat_inv * Phi_pp * Z_mat) - (Phi_up*Z_mat - Z_mat_inv*Phi_pu) );
Phi_22_hat_true = 0.5*( (Phi_uu + Z_mat_inv * Phi_pp * Z_mat) + (Phi_up*Z_mat + Z_mat_inv*Phi_pu) );

norm(Phi_11_hat_true, 1)
norm(Phi_11_hat_true, 2)
norm(Phi_11_hat_true, inf)



Delta_11 = Phi_11_hat_true - Phi_hat_11;
Delta_12 = Phi_12_hat_true - Phi_hat_12;
Delta_21 = Phi_21_hat_true - Phi_hat_21;
Delta_22 = Phi_22_hat_true - Phi_hat_22;

norm(Delta_11, inf)
norm(Delta_12, inf)
norm(Delta_21, inf)
norm(Delta_22, inf)
% 
% full(Phi_12_char_true)
% full(Phi_char_12)

%full(Phi_11_char_true - Phi_char_11)