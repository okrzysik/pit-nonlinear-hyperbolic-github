clc
clear

k = 2;

[zeta_left, zeta_right] = reconstruction_weights(k);
[d_left,    d_right]    = linear_weights(k);

% -1 * zeta_right.*d_right * factorial(2*k) / (factorial(k+1) * factorial(k) * factorial(k-1))
% 
% +1 * zeta_left.*d_left * factorial(2*k) / (factorial(k+1) * factorial(k) * factorial(k-1))

zeta_right.*d_right
zeta_left.*d_left

function [zeta_left, zeta_right] = reconstruction_weights(k)

    zeta_left = [];
    zeta_right = [];

    for ell = 0:k-1
        zeta_left  = [zeta_left;  (-1)^(ell)*factorial(ell)*factorial(k-ell)      ];
        zeta_right = [zeta_right; (-1)^(ell+1)*factorial(ell+1)*factorial(k-ell-1)];
    end

end


function [d_left, d_right] = linear_weights(k)
    if k == 1
        d_right = 1;
    elseif k == 2
        d_right = [2/3; 1/3];
    elseif k == 3
        d_right = [3/10; 3/5; 1/10];
    end

    d_left = flipud(d_right);
end

% for k = 1:9
% %factorial(k) * factorial(k-1) / factorial(2*k)
% 
% factorial(2*k) / (factorial(k+1) * factorial(k) * factorial(k-1))
% 
% end