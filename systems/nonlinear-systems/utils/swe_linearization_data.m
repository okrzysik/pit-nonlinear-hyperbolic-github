% The class implements data structures for storing linearized information
% from the SWE. This is the data needed to solve the space-time 
% problem via global linearization.
classdef swe_linearization_data < linearization_data
    
    properties
       h_pos;
       h_neg;
       u_pos;
       u_neg;
    end
    
    methods 
        function obj = swe_linearization_data(nt_)
            % Call parent's constructor
            obj = obj@linearization_data(nt_);
            
            obj.h_pos = cell(obj.nt, 1);
            obj.h_neg = cell(obj.nt, 1);
            obj.u_pos = cell(obj.nt, 1);
            obj.u_neg = cell(obj.nt, 1);
        end
    end
    
end