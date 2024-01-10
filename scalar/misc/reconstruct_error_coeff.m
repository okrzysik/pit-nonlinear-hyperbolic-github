% The zeta coefficients that weight the derivative in the error estimate
% for reconstruction at interfaces.
%
% There's also a check that the factorial-based expressions for the
% coefficients are consistent with the original formula I obtained (they
% are).

clc
k = 2; % The polynomial recovers cell averages for k cells/it's primitive 
% interpolates the primitive at k+1 cell interfaces

fprintf('Reconstruction error coefficients, k=%d\n', k)

fprintf('LHS interface\n')
for ell = 0:k-1
    
    zeta_left = (-1)^(ell)*factorial(ell)*factorial(k-ell);
    fprintf('ell=%d, zeta = %d\n', ell, zeta_left);
    
%     t = 1;
%     for j = -ell-1:k-ell-1
%        
%         if j ~= -1
%             t = t * (j+1);
%         end
%     end
%     [t, zeta] 
end

fprintf('RHS interface\n')
for ell = 0:k-1
    zeta_right = (-1)^(ell+1)*factorial(ell+1)*factorial(k-ell-1);
    fprintf('ell=%d, zeta = %d\n', ell, zeta_right);
    
%     t = 1;
%     for j = -ell-1:k-ell-1
%        
%         if j ~= 0
%             t = t * j;
%         end
%     end
%     [t, zeta]
end


