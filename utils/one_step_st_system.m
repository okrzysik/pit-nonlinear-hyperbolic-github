% Implements abstract class for space-time system A(q) = b which represents
% the system of equations,
%   q_{1}   = b_{1}.
%   q_{n+1} = Phi(q_{n}) + b_{n+1}, n = 1, ..., n_t - 1.
%
% Note: There are nt time points, and nt equations. Each of q and b has nt
% blocks
%
% This class can perform several key functions on the system including
% residual computations and CF-based relaxations.

classdef one_step_st_system < handle
   
    properties
        step             % Handle for the step function, Phi. Must conform to q1 = step(t0_idx, q0)
        step_dimension   % Int. Length of vectors that go in and out of Phi. This is fixed independent of the time-step size.
        nt               % Int. Number of points on time grid.

        CF_splitting_pa % Parameter struct with information about CF-based splitting and associated relaxation.
                        % Requires the fields:
                        %    cf.           INT.    The number determining which the CF-splitting, needs to be bigger than 1. Every cf-th point, starting at 1, is a C-point.
                        %    relax_scheme. STRING. CF-based relaxation scheme, e.g., 'F', 'FCF', 'FCF', etc.
    end
    
   
    methods
        
    %% Constructor
    function obj = one_step_st_system(step_, step_dimension_, nt_, CF_splitting_pa_)
        
        obj.step           = step_;
        obj.step_dimension = step_dimension_;
        obj.nt             = nt_;
        
        if nargin == 3
            obj.CF_splitting_pa = [];
        else
            obj.CF_splitting_pa = CF_splitting_pa_;
        end
    end
        
    %% Solve the system by forward substitution, i.e., time-stepping
    % to solve the system A*q = b.
    function q = forward_solve(obj, b)
        b = reshape(b, [obj.step_dimension, obj.nt]);
        q = zeros(obj.step_dimension, obj.nt); 
        q(:, 1) = b(:, 1); % Initial condition.
        for t0idx = 1:obj.nt-1
            q(:, t0idx+1) = obj.step(t0idx, q(:, t0idx)) + b(:, t0idx + 1);
        end
        q = q(:);
    end
    
    %% Compute the action of A on some vector c: b = A(v)
    function c = system_action(obj, v)
        v = reshape(v, [obj.step_dimension, obj.nt]);
        c = zeros(size(v));
        
        c(:, 1) = v(:, 1);
        for tidx = 1:obj.nt-1
            c(:, tidx+1) = v(:, tidx+1) - obj.step(tidx, v(:, tidx)); % (A(q))_{n+1} = q_{n+1} - Phi(q_{n}).
        end
        
        c = c(:);
    end
    
    %% Space-time residual r = b - A(q).
    function r = residual(obj, q, b)
        b = reshape(b, [obj.step_dimension, obj.nt]);
        q = reshape(q, [obj.step_dimension, obj.nt]);
        r = zeros(size(q));
        r(:, 1) = b(:, 1) - q(:, 1); % Initial condition.
        for t0idx = 1:obj.nt-1
            r(:, t0idx+1) = b(:, t0idx+1) - ( q(:, t0idx+1) - obj.step(t0idx, q(:, t0idx)) ); % (A(q))_{n+1} = q_{n+1} - Phi(q_{n}).
        end
        
        r = r(:);
    end
    
    
    %% Relaxation for space-time system. r is the residual, and is only 
    % computed if requested. The residual is r = b - A(q).
    function [q, r] = relaxation(obj, q, b)
        
        btemp = reshape(b, [obj.step_dimension, obj.nt]);
        qtemp = reshape(q, [obj.step_dimension, obj.nt]);

        for relax_step_idx = 1:numel(obj.CF_splitting_pa.relax_scheme)
        
            cf = obj.CF_splitting_pa.cf;
            assert(cf ~= 1, 'Ambiguous what CF splitting is for cf = 1.')
            
            relax_step = obj.CF_splitting_pa.relax_scheme(relax_step_idx);
            
            if strcmp(relax_step, 'F')
                qtemp = F_relax(obj, qtemp, btemp, cf);
                
            elseif strcmp(relax_step, 'C')
                qtemp = C_relax(obj, qtemp, btemp, cf);
                
            else
                error('Relax step = %s not recognised... Must be C or F.', relax_step)
            end
        end

        % Compute the residual if it was requested. This is based on what
        % the last relaxation step was. C-relax leaves F-points alone, and
        % F-relax leaves C-points alone. 
        % We can use this to compute the full residual as below.  
        if nargout == 2
            if exist('relax_step', 'var')
                if strcmp(relax_step, 'F')
                    r = C_relax(obj, qtemp, btemp, cf) - qtemp;
                elseif strcmp(relax_step, 'C')
                    % This residual computation is a little bit wasteful
                    % since it doesn't exploit the fact that residuals at
                    % C-points are zero, but typically we'd have a
                    % coarsening factor >> 1, so that the wasted
                    % computation here is not significant.
                    r = residual(obj, qtemp(:), b(:));
                end
                r = r(:);
            % No relaxation scheme, so compute full residual vector
            else
                r = residual(obj, q(:), b(:));
            end
        end
        
        q = qtemp(:); 
    end

    function q = F_relax(obj, q, b, cf)
        for CF_interval = 1:ceil(obj.nt/cf)
            for F_pt_idx = 1:cf-1
                t0_idx = (CF_interval-1)*cf + F_pt_idx;
                if t0_idx > obj.nt - 1; break; end
                q(:, t0_idx+1) = b(:, t0_idx+1) + obj.step(t0_idx, q(:, t0_idx));
            end    
        end
    end

    function q = C_relax(obj, q, b, cf)
        q(:, 1) = b(:, 1); % Initial condition.
        for CF_interval = 1:ceil(obj.nt/cf)
            t0_idx = CF_interval*cf;
            if t0_idx > obj.nt-1; break; end
            q(:, t0_idx+1) = b(:, t0_idx+1) + obj.step(t0_idx, q(:, t0_idx));
        end
    end
    %% End of relaxation
    
    %% Update all-point values of the vector q by applying the transform T to them.
    function Tq = all_point_transformation(obj, q, T)
        
        qtemp = reshape(q, [obj.step_dimension, obj.nt]);
        Tqtemp = qtemp; % Initialize as same vector since only some entries are over-written.
        
        for tidx = 1:obj.nt
            Tqtemp(:, tidx) = T(tidx, qtemp(:, tidx));
        end
        
        Tq = Tqtemp(:);
    end
    
    %% Update C-point values of the vector q by applying the transform T to them.
    % T == T(t0_idx, q(t0_idx))
    function Tq = C_point_transformation(obj, q, T)
        
        cf = obj.CF_splitting_pa.cf;
        assert(cf ~= 1, 'Ambiguous what CF splitting is for cf = 1.')
        qtemp = reshape(q, [obj.step_dimension, obj.nt]);
        Tqtemp = qtemp; % Initialize as same vector since only some entries are over-written.
        
        % The first point is a C-point
        Tqtemp(:, 1) = T(1, qtemp(:, 1));
        
        for CF_interval = 1:ceil(obj.nt/cf)
            t0_idx = CF_interval*cf;
            if t0_idx > obj.nt-1; break; end
            Tqtemp(:, t0_idx+1) = T(t0_idx+1, qtemp(:, t0_idx+1));
        end
        
        Tq = Tqtemp(:);
    end
    
    
    %% Update C-point values of the vector q by applying the transform T to them.
    % T == T(t0_idx, q(t0_idx))
    function Tq = F_point_transformation(obj, q, T)
        
        cf = obj.CF_splitting_pa.cf;
        assert(cf ~= 1, 'Ambiguous what CF splitting is for cf = 1.')
        qtemp  = reshape(q, [obj.step_dimension, obj.nt]);
        Tqtemp = qtemp; % Initialize as same vector since only some entries are over-written.
        
        for CF_interval = 1:ceil(obj.nt/cf)
            for F_pt_idx = 1:cf-1
                t0_idx = (CF_interval-1)*cf + F_pt_idx;
                if t0_idx > obj.nt - 1; break; end
                Tqtemp(:, t0_idx+1) = T(t0_idx+1, qtemp(:, t0_idx+1));
            end
        end
        
        Tq = Tqtemp(:);
    end
        
        
    end
    
end