% Acoustics material parameters c0 and Z0.
%   
% The parameters will be returned evaluated at the x that is parsed. As a
% minimum, the struct "pde_pa" must contain the field "mat_param_id"
%
% The layered medium examples must also contain the field
% "mat_param_num_layer" and can contain the optional field 
% "mat_param_homogenize" which homegenizes the medium to have constant
% sound speed and impedance (by default these problems are not homogenized) 
%
%
% Examples 1,2,3 come from Section 5 of 
% "A WAVE PROPAGATION METHOD FOR CONSERVATION LAWS AND BALANCE LAWS WITH 
% SPATIALLY VARYING FLUX FUNCTIONS" by Bale et al. (2002)
%
% Typos in Bale et al. examples (as far as I can tell, anyway):
%   Ex. 1, Fig 4: Solution in middle pannel is shown at t = 0.4 and NOT t = 0.3 as claimed
%   Ex. 2, Fig 5: Solution in middle pannel is shown at t = 0.4 and NOT t = 0.35 as claimed
%
%   Ex. 3, As written in Bale et al.'s (5.8) is wrong... Just look at the 
%   left plot in the corresponding fig. The intervals should be reversed.
function [c, Z] = acoustics_material_parameters(x, pde_pa)

    if pde_pa.mat_param_id == 0
        c = ones(size(x));
        Z = ones(size(x));

    elseif pde_pa.mat_param_id == 1
        c = 1 + 0.5 * sin(10*pi*x);
        Z = ones(size(x));

    elseif pde_pa.mat_param_id == 2
        c = 1 + 0.5 * sin(10*pi*x);
        Z = 1 + 0.25 * cos(10*pi*x);
        %Z = 1 + 0.75 * cos(10*pi*x);

    elseif pde_pa.mat_param_id == 3
        
        I = find(x > 0.35 & x < 0.65);
        c = 0.6*ones(size(x));
        Z = 6*ones(size(x));
        c(I) = 2;
        Z(I) = 2;
    
    elseif pde_pa.mat_param_id == 4
        I = find(x > 0.35 & x < 0.65);
        c = ones(size(x));
        Z = ones(size(x));
        c(I) = 0.5;
        Z(I) = 2;

    % Periodically repeating layered medium. See LeVeque (2004), example
    % 9.4.
    elseif pde_pa.mat_param_id == 5
        
        
        num_layers       = pde_pa.mat_param_num_layers;
        points_per_layer = numel(x)/num_layers;
        
        assert(floor(points_per_layer) == points_per_layer, ...
            'Set up requires integer # grid points per layer: %d (# points) not divisible by %d (# layers)', numel(x), num_layers)

        
        if isfield(pde_pa, 'mat_param_homogenize')
            homogenize = pde_pa.mat_param_homogenize;
        else
            homogenize = ~true;
        end
        
        epsilon = 1;
        rho_odd_layer  = 1;
        K_odd_layer    = 1;
        rho_even_layer = 1 + epsilon;
        K_even_layer   = 1 + epsilon;
        
        if ~homogenize
            for layer_idx = 1:num_layers
                for j = 1:points_per_layer
                    point_idx = (layer_idx-1)*points_per_layer + j;
                    if mod(layer_idx, 2) == 1
                        rho(point_idx) = rho_odd_layer;
                        K(point_idx)   = K_odd_layer;
                    else
                        rho(point_idx) = rho_even_layer;
                        K(point_idx)   = K_even_layer;
                    end
                end
            end
        elseif homogenize
           rho = (rho_even_layer + rho_odd_layer)/2            * ones(size(x));
           K   = (1/2 * (1/K_even_layer + 1/K_odd_layer) )^-1  * ones(size(x));
        end
        
        c = sqrt(K ./ rho);
        Z = sqrt(K .* rho);
    
        c = c(:);
        Z = Z(:);
        
    % Randomly varying layered medium. 
    % See: "High-resolution finite-volume methods for acoustic
    % waves in periodic and random media" by Fogarty and LeVeque
    % (1999).
    % From F & L:         
    % Coefficients vary from 0.2 up to 1.8. rho and K are taken to be 
    % piecewise constant. We use a layered 
    % medium with 2000 layers in which values of rho and K are randomly 
    % chosen from uniform distributions between 0.2 and 1.8.
    elseif pde_pa.mat_param_id == 6
            
        num_layers       = pde_pa.mat_param_num_layers;
        points_per_layer = numel(x)/num_layers;
        
        assert(floor(points_per_layer) == points_per_layer, ...
            'Set up requires integer # grid points per layer: %d (# points) not divisible by %d (# layers)', numel(x), num_layers)

        if isfield(pde_pa, 'mat_param_homogenize')
            homogenize = pde_pa.mat_param_homogenize;
        else
            homogenize = ~true;
        end
        
        rho_background = 0.2;
        K_background   = 0.2;
        rng(51); % Seed RNG here so that the same random values are 
        % always chosen independent of when this function is called.
        rho_noise = 1.6*rand(num_layers, 1);
        K_noise   = 1.6*rand(num_layers, 1);

        if ~homogenize
            for layer_idx = 1:num_layers
                for j = 1:points_per_layer
                    point_idx      = (layer_idx-1)*points_per_layer + j;
                    rho(point_idx) = rho_background + rho_noise(layer_idx);
                    K(point_idx)   = K_background + K_noise(layer_idx);
                end
            end
        elseif homogenize
            rho_mean    = mean(rho_background + rho_noise);
            rho(:)      = rho_mean;

            K_harm_mean = (1/numel(x)*( sum(1./(K_background + K_noise)) ))^-1;
            K(:)        = K_harm_mean;
        end
        
        c = sqrt(K ./ rho);
        Z = sqrt(K .* rho);
    
        c = c(:);
        Z = Z(:);
        
    else
        error('Material parameters not recognised: pde_pa.mat_param_id = %d', pde_pa.mat_param_id)
    end
end