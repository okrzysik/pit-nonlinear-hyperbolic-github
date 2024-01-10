% Compute an itegral via high-order Gauss--Legendre quadrature.
function ubar = gauss_legendre_quad(u, xl, xr, n)

    if n ~= 7
        error('Only 7-point Gauss--Legendre quadrature implemented')
    end

    % Integration nodes in [-1, 1]
    unit_interval_nodes = [ -0.949107912342759;
                            -0.741531185599394;
                            -0.405845151377397;
                             0.000000000000000;
                             0.405845151377397;
                             0.741531185599395;
                             0.949107912342758];
    % Map nodes from [-1, 1] into [xl, xr]       
    scaled_nodes = 0.5*(xr - xl)*unit_interval_nodes + 0.5*(xr + xl);
    
    weights = [ 0.129484966168870;
                0.279705391489277;
                0.381830050505119;
                0.417959183673470;
                0.381830050505119;
                0.279705391489278;
                0.129484966168869];
    
    % Compute sum via an inner product and multiply by original domain. 
    ubar = u(scaled_nodes)' * weights * 0.5 * (xr - xl);
end

% I got these from weights and nodes from an existing code 
% /Users/oliverkrzysik/Library/CloudStorage/OneDrive-UniversityofWaterloo/projects/nonlinear_st/c/misc/gauss_quad_rules.c
% and see also the code 
% /Users/oliverkrzysik/Library/CloudStorage/OneDrive-UniversityofWaterloo/projects/nonlinear_st/python/miscellaneous/numerical_cell_averages.py
% which shows how to compute evaluate the quadrature rule.