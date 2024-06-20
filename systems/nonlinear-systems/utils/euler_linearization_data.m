% The class implements data structures for storing linearized information
% from the Euler equations. This is the data needed to solve the space-time 
% problem via global linearization.
classdef euler_linearization_data < linearization_data 
    
    properties
       rho_pos;
       rho_neg;
       u_pos;
       u_neg;
       E_pos;
       E_neg;
    end
    
    methods 
        function obj = euler_linearization_data(nt_)
            % Call parent's constructor
            obj = obj@linearization_data(nt_);
            
            obj.rho_pos = cell(obj.nt, 1);
            obj.rho_neg = cell(obj.nt, 1);
            obj.u_pos   = cell(obj.nt, 1);
            obj.u_neg   = cell(obj.nt, 1);
            obj.E_pos   = cell(obj.nt, 1);
            obj.E_neg   = cell(obj.nt, 1);
        end
    end
    
end