% Data structure that stores linearization data required for the space-time
% solution of a linear conservation law.
%
classdef scalar_linearization_data < linearization_data
   
    properties
        cell_averages           
        reconstruction_weights 
        dissipation            
        wave_speed
    end
    
    methods 
        
        function obj = scalar_linearization_data(nt_)
            % Call parent's constructor
            obj = obj@linearization_data(nt_);
            
            obj.cell_averages          = cell(nt_, 1);
            obj.reconstruction_weights = cell(nt_, 1);
            obj.dissipation            = cell(nt_, 1);
            obj.wave_speed             = cell(nt_, 1);
        end
    end
    
end