function C = stencil2matrix(stencil, diags, dimC, circulant_flag)
%CIRCULANT assemble sparse banded circulant matrix:
%
%INPUT: 
%   stencil             :   Array of doubles. Stencil of matrix C. Ordered from lowest non-zero
%                               diagonal up to the highest non-zero diagonal.
%   diags               :   Array of ints. Diagonal indices of the elements in stencil.
%                               Main diagonal has index 0, subdiagonals have
%                               index < 0, superdiagonals have index > 0.
%   n                   :   Int. Dimension of (square) matrix
%   circulant_flag      :   INT. OPTIONAL. Default is 1 for circulant
%                               matrix. Pass as 0 for non-circulant Toeplitz matrix.
%
%OUTPUT:
%   C           :   Sparse banded circulant matrix
%
%NOTES:
%   Too hard to write out the banded circulant matrix above. You get the
%   point.
%
%
%TODO:
% There can be some inconsistency here I think, where I can overspecify a
% circulant matrix in terms of specifying its first and last diagonals (this could even happen implicitly I think when I compute things below)...
% That's dumb... I should think more about this... There needs to be some
% limitation on how close I can specify a diagonal close to the boundary
% relative to the other boundary...

%%
% Check for non-circulant flag. Build circulant matrix if nothing is passed
if nargin < 4
    circulant_flag = 1;
elseif isempty(circulant_flag)
    circulant_flag = 1;
elseif circulant_flag ~= 1 && circulant_flag ~= 0
    error('Pass 0 for non-circulant Toeplitz. Pass 1 for circulant Toeplitz')
end

% Just check that diags has only integers
if sum( sum(abs(diags - floor(diags))) ) > 0
    error('diags must be an array of integers...')
end

% Check information passed is consistent
if numel(diags) ~= numel(stencil)
    error('diags and stencil must have same number of entries')
end

if min(diags) < -dimC+1 || max(diags) > dimC-1
    
    
    
    error('diagonal indices diags must not exceed the dimension of the matrix')
end

% Write as column vectors in case they're not already.
stencil = stencil(:); 
diags = diags(:);



%% Build banded non-circulant Toeplitz matrix in sparse format.
if circulant_flag == 0
    % Sort diagonal indices, maybe spdiags appreciates this?
    [diags, idx] = sort(diags);
    stencil = stencil(idx);
    
    e = ones(dimC, 1);
    C = spdiags(e*stencil', diags, dimC, dimC);

    
%% Build banded circulant matrix in sparse format.
elseif circulant_flag == 1

    % Sort diagonal indices, maybe spdiags appreciates this?
    [diags, idx] = sort(diags);
    stencil = stencil(idx);

    % Identify indices of sub and super diagonals
    subdiags_idx   = find( diags < 0 );
    superdiags_idx = find( diags > 0 );

    % Get sub and super diagonal indices
    subdiags   = diags(subdiags_idx);
    superdiags = diags(superdiags_idx);
    
    % Get entries of sup and super diagonals
    substencil   = stencil(subdiags_idx);
    superstencil = stencil(superdiags_idx);
    
    % Encorporate these extra terms into the stencil and diags arrays, wrapping them around by n
    diags = [superdiags-dimC; diags; dimC+subdiags];
    stencil = [superstencil; stencil; substencil];
   
    e = ones(dimC, 1);
    C = spdiags(e*stencil', diags, dimC, dimC);
end

end