% Decompose the arrival information of the FV cell.
function [integral_translation, epsilon] = SL_decompose_FV_cell_evolution(arrive, depart, h)
        
        % The mesh-normalized distance from the arrival point to the
        % departure point. The ceiling of this number is the number of
        % cells away the east boundary of the departure point cell is from
        % the east boundary that is the arrival point. The remainder 
        % is the mesh-normalized distance that the
        % departure point is from the east boundary of the cell that it's in.
        mesh_normalized_dist = (depart - arrive)/h;

        % Number of cells that departure point is translated from arrival 
        % point (as measured from the east boundary of the cell containing the departure point). 
        integral_translation = ceil(mesh_normalized_dist);
        
        % Distance from departure point to the cell boundary to its east.
        epsilon = integral_translation - mesh_normalized_dist;
end