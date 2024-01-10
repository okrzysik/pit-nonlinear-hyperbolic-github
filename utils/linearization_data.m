% Data structure for storing linearization data. Since this is a handle
% subclass it will effectively act as though it's a pointer.
%
% Different linearization data is required for different
% PDEs/discretizations so this class is intended to act as a parent class.
classdef linearization_data < handle
   
    properties
        nt        = []; % Number of time points
        t0_idx    = []; % Index where linearization data will be accessed/written
        stage_idx = []; % Stage index for multi-stage RK methods
    end
    
    methods 
       function obj = linearization_data(nt_)
            obj.nt = nt_;
       end
       
       function set.t0_idx(obj, t0_idx_)
           obj.t0_idx = t0_idx_;
       end
       
       function set.stage_idx(obj, stage_idx_)
           obj.stage_idx = stage_idx_;
       end
    end
    
end