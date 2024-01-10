function [u_vector, rnorms, MGRIT_object, enorms] = mgrit(u_vector, g_vector, step, MGRIT_object, solver_params)
%
% Use MGRIT to solve the linear system: A*u = g. 
%
% u_vector: initial iterate vector
%
% g_vector: right-hand side vector
%
% step: HANDLE. Applies the time-stepping operator. Must be setup as
%       [u1, MGRIT_object] = step(u0, step_status, MGRIT_object), 
%   where:
%       u0: The vector to be stepped from
%       step_status: STRUCT holding t0_idx, t0, dt, and level.
%       MGRIT_object: STRUCT. See below.
%       u1: The result of applying the time-stepping operator to u0.
%
% MGRIT_object: STRUCT. This object is parsed in and out of every call to 
% STEP and in and out of the call to MGRIT. Solver-specific information is 
% stored here, including information about the hierarchy. The user is free 
% to put information in this, but they shouldn't overwrite any data in it
% that they didn't put there.
%
% solver_params: STRUCT. Holds all of the required solver parameters. See
%   code below for what's required.
%
%
% NOTES:
% - Linear MGRIT is implemented, and NOT FAS-MGRIT. Linear MGRIT is more
% efficient than FAS-MGRIT when applied to linear problems since per
% two-grid cycle, FAS requires the coarse-grid operator be applied twice. 
%
% - Levels are indexed such that the finest level has index 1, the next
% coarsest level has index 2, etc.
%
% - A two-level MGRIT iteration is:
% 1. u <-- relax(u) Pre-relaxation, this is F, FCF, FCFCF, etc.
% 2. r_c <-- R(g - A*u) Compute the coarse-grid residual. Since pre-relax ends
%       with an F-relaxation, this just means computing the residual at C-points
% 3. u <-- u + inv(A_c)*r_c Coarse-grid correction
% 4. if terminating, possibly do an F-relaxation
%
% - The residual norms that are reported back are computed based on those in
% step 2. of the above algorithm. Thus, the initial residual norm is not
% really the residual norm of the initial iterate. Furthermore, the final
% residual norm that is reported is not actually that of the final iterate
% because the iterate has since had a coarse-grid correction applied to it,
% and possibly an F-relaxation. 


%% Check and adjust inputs as required

% Required solver parameters.
sp_required = {...
    'final_F_relax', ...    % BOOLEAN. Apply final F-relaxation at end of solve phase.
    'cf', ...               % INT or level-dependent HANDLE. Coarsening factor. 
    'pre_relax', ...        % STRING or level-dependent HANDLE. Level-dependent relaxation. 
    'res_halt_tol', ...     % DOUBLE. Absolute halting tolerance for residual norm.
    'maxlevels', ...        % INT. Maximum number of levels in multigrid hierarchy.
    'maxiter', ...          % INT. Maxmimum number of MGRIT iterations. 
    'min_coarse_nt', ...    % INT. Minimum number of points allowed on the coarsest level.
};

% Optional solver parameters and their default values.
sp_optional = {...
    {'res_reduction'; true}; ...          % Halt based on the residual norm being reduced by a certain amount rather than its absolute value.
    {'warm_restart';  false}; ...         % If a hierarchy already exists in MGRIT_object and it is to be used here
    {'verbose';       true}; ...          % Spit out MGRIT messages or not
    {'verbose_assembly'; ~true}; ...      % Spit out MGRIT messages or not about assembling hierarchy
    {'res_norm_handle'; @(r) norm(r, 2)}; % Norm used for monitoring residual (and error)
};

% Basic check that all mandatory fields exisit in solver_params
for spr_idx = 1:numel(sp_required)
    if ~isfield(solver_params, sp_required{spr_idx})
       error('solver_params requires field: %s\n', sp_required{spr_idx}) 
    end
end

% Add optional values into solver params if they didn't exist already
for spo_idx = 1:numel(sp_optional)
    if ~isfield(solver_params, sp_optional{spo_idx}{1})
       solver_params.(sp_optional{spo_idx}{1}) = sp_optional{spo_idx}{2};
    end
end

% If level-dependent options were parsed as single objects, then create 
%level-dependent functions extrapolating those values to any level
if ~isa(solver_params.pre_relax, 'function_handle')
   pre_relax1 = solver_params.pre_relax;
   solver_params.pre_relax = @(level) pre_relax1;
end

if ~isa(solver_params.cf, 'function_handle')
   cf1 = solver_params.cf;
   solver_params.cf = @(level) cf1;
end


if nargout == 4
    compute_errors = true;
else
    compute_errors = false;
end


%% Initialize things required before solving.

% Assemble hierarchy
if isfield(MGRIT_object, 'hierarchy') && ~isempty(MGRIT_object.hierarchy) && solver_params.warm_restart
    
else
     MGRIT_object.hierarchy = assemble_hierarchy(MGRIT_object, solver_params);
end

level  = 1;
u_cell = vector2cell(u_vector, level, MGRIT_object);
g_cell = vector2cell(g_vector, level, MGRIT_object);

%% If errors requested then compute the exact solution
if compute_errors 
   u_exact_cell = time_stepping_solve(g_cell, 1, step, MGRIT_object); 
   e_cell       = cellfun(@minus, u_exact_cell, u_cell, 'Un', 0); % https://www.mathworks.com/matlabcentral/answers/8793-subtracting-two-cell-arrays-yielding-a-third-cell-array
   enorms       = solver_params.res_norm_handle(cell2mat(e_cell));
end

% Print progress to command line
% Execute MGRIT V cycles
if solver_params.verbose
    fprintf('--------------------------------------------\n')
    fprintf(' --------- Executing MGRIT cycles --------- \n')
    fprintf('--------------------------------------------\n')
    fprintf('| iter |     |r|    |  |r|/|r0|  |   conv.   |\n')
    fprintf('----------------------------------------------\n')
end


%% Apply V-cycles
iters = 1;
rc_norm = Inf; % Set to a dummy value, ensuring we do at least 1 iteration.

% Upack residual halting tolerance.
rtol = solver_params.res_halt_tol;

while rc_norm > rtol && iters <= solver_params.maxiter

    level  = 1;
    MGRIT_object.iter = iters;
    [u_cell, MGRIT_object, rc_norm] = vcycle_linear(u_cell, g_cell, level, step, MGRIT_object, solver_params);
    rnorms(iters) = rc_norm;
    
    % If we're going to monitor residual reduction, then scale halting 
    % tolerance by initial residual norm.
    if iters == 1 && solver_params.res_reduction
        rtol = rtol * rc_norm;
    end
    
    if compute_errors 
       e_cell = cellfun(@minus, u_exact_cell, u_cell, 'Un', 0); % https://www.mathworks.com/matlabcentral/answers/8793-subtracting-two-cell-arrays-yielding-a-third-cell-array
       enorm  = solver_params.res_norm_handle(cell2mat(e_cell));
       enorms(iters) = enorm;
    end

    % Print residual monitoring to command line
    if solver_params.verbose
        if iters == 1
            fprintf('|  0   |  %.2e  |  %.2e  |    ---    |\n', rnorms(1), 1);
        else
            if iters > 10; dum = ''; else; dum = ' '; end
            fprintf('|  %d%s  |  %.2e  |  %.2e  |  %.2e |\n', iters-1, dum, rnorms(iters), rnorms(iters)/rnorms(1), rnorms(iters)/rnorms(iters-1));
        end
    end
    
    iters  = iters + 1; % Increment iteration counter
end




%% Before exiting solve, apply final F-relaxation if requested by user.
if solver_params.final_F_relax
    [u_cell, MGRIT_object] = relaxation('F', u_cell, g_cell, level, step, MGRIT_object);
end

% Convert the solution back into a vector.
u_vector = cell2mat(u_cell);

end


%% Sub functions

%% Create a time-grid hierarchy based on t.
function hierarchy = assemble_hierarchy(MGRIT_object, solver_params)
    hierarchy = struct([]);

    level = 1;

    % Build time grid hierarchy.
    while level <= solver_params.maxlevels

        m = solver_params.cf(level);
        
        % Define coarsest level to use a coarsening factor of 1.
        if level == solver_params.maxlevels
            m = 1;
        end

        % On coarse levels, grid comes from coarsening parent grid.
        if level > 1
            t = hierarchy(level-1).t(hierarchy(level-1).cpoint_indptr(1:end-1));
        else
            t = MGRIT_object.t;
            assert(numel(t) > solver_params.min_coarse_nt, "fine grid has more points that minimum allowable number");
        end
            
        nt = numel(t);

        % Drop out of the loop if the current grid has too few points.
        if nt < solver_params.min_coarse_nt
            break
        end

        hierarchy(level).t = t;

        % Number of time points on current level.
        hierarchy(level).nt = numel(t);

%         % It doesn't make sense to set up coarse grid if already on the coarsest grid
%         if level ~= solver_params.maxlevels

        % Get relaxation scheme.
        hierarchy(level).pre_relax  = solver_params.pre_relax(level);

        % Compute indices of C-points on current grid by coarsening every mth point. 
        cpoint_inds = [1:m:numel(t)].';

        % If the coarsened grid has too few points then don't actually coarsen.
        if numel(cpoint_inds) < solver_params.min_coarse_nt 
            break
        end

        hierarchy(level).cpoint_indptr = [cpoint_inds; nt+1];

        % Store coarsening factor that defines child grid.
        hierarchy(level).cf = m;
        %end

        % Go to next level in hierarchy.
        level = level + 1;
    end
    

    % Print out hierarchy statistics
    if solver_params.verbose_assembly || solver_params.verbose
       fprintf('\n---------------------------------------------\n');
       fprintf('| ------------- MGRIT hierarchy -------------- |\n');
       fprintf('-----------------------------------------------\n');
       fprintf('| Level |    nt    | total DOFs |  m  |\n');
       fprintf('----------------------------------------------\n');
       for level = 1:numel(hierarchy)
           fprintf('|   %d   | %.2e |  %.2e |  %d  |\n', ...
               level, hierarchy(level).nt, ...
               hierarchy(level).nt * MGRIT_object.block_size, hierarchy(level).cf); 
       end
       fprintf('----------------------------------------------\n\n');
    end

% End build_hierarchy    
end

%% Map a vector into a time-dependent cell array
function u_cell = vector2cell(u_vector, level, MGRIT_object)

    nt         = MGRIT_object.hierarchy(level).nt;
    u_cell     = cell(nt, 1);

    % Unpack vectors into cells.
    for n = 1:nt
        start = (n-1)*MGRIT_object.block_size+1;
        stop  = n*MGRIT_object.block_size;
        u_cell{n} = u_vector(start:stop);
    end
end
% End of vector2cell



%% Apply a V-cycle. 
% The norm of the coarse-grid residual is returned 
function [u_cell, MGRIT_object, rc_norm] = vcycle_linear(u_cell, g_cell, level, step, MGRIT_object, solver_params)

    coarsest_level = numel(MGRIT_object.hierarchy);

    % Recursively solve coarse-grid problem.
    if level < coarsest_level

        % Unpack some components to save typing.
        nt_coarse          = MGRIT_object.hierarchy(level+1).nt; % Number of time points on coarse level
        cpoint_inds_global = MGRIT_object.hierarchy(level).cpoint_indptr; % Global indices of C-points on the current grid.

        % Pre-relaxation
        pre_relax  = MGRIT_object.hierarchy(level).pre_relax;  % String
        %assert(strcmp(pre_relax(1), 'F'), 'pre-relaxation must begin with an F-relax!')
        
        for i = 1:numel(pre_relax)
            [u_cell, MGRIT_object] = relaxation(pre_relax(i), u_cell, g_cell, level, step, MGRIT_object);
        end
        
        assert(strcmp(pre_relax(i), 'F'), 'pre-relaxation must end with an F-relax so that restriction is ideal!')
        
        % Compute coarse-grid residual
        [rc_cell, MGRIT_object] = Cpoint_residual(u_cell, g_cell, level, step, MGRIT_object);
        % Compute norm of this coarse-grid residual if on the finest level
        % (for residual monitoring)
        if level == 1
            rc_norm = solver_params.res_norm_handle(cell2mat(rc_cell));
        end

        % Recurse to coarse grid to compute C-point error, use 0 initial guess.
        ec_cell = cell(nt_coarse, 1);
        for idx = 1:nt_coarse
            ec_cell{idx} = zeros(MGRIT_object.block_size, 1);
        end
        [ec_cell, MGRIT_object] = vcycle_linear(ec_cell, rc_cell, level+1, step, MGRIT_object);

        % Correct solution at C-points using injection interpolation.
        for cpoint_idx_local = 1:nt_coarse
            u_cell{cpoint_inds_global(cpoint_idx_local)} = ...
                u_cell{cpoint_inds_global(cpoint_idx_local)} + ec_cell{cpoint_idx_local};
        end
        
        % Post F-relaxation.
        % If on level 1, no need to apply post F-relaxation because the
        % pre-relaxation on the next iteration will apply an F-relaxation
        % (the execption being if there's no next iteration, but this is
        % covered outside of the V-cycle code).
        if level == 1
            
        else
            [u_cell, MGRIT_object] = relaxation('F', u_cell, g_cell, level, step, MGRIT_object);
        end

    % On coarsest level, use alternative solver here.
    elseif level == coarsest_level
        [u_cell, MGRIT_object] = time_stepping_solve(g_cell, level, step, MGRIT_object);
        rc_norm = inf; % There is no C-point residual in this case.
    end

end
% End of vcycle_linear


%%
function [u_out_cell, MGRIT_object] = relaxation(scheme, u_in_cell, g_cell, level, step, MGRIT_object)

    step_status.level = level;
    step_status.parent_function = 'relaxation';

    t                  = MGRIT_object.hierarchy(level).t;
    nt_coarse          = MGRIT_object.hierarchy(level+1).nt;          % Number of coarse-grid time points.
    cpoint_global_inds = MGRIT_object.hierarchy(level).cpoint_indptr; % Global indices of C-points on the current grid.

    % Initialize output array.
    u_out_cell = u_in_cell;
    
    MGRIT_object.u_current = u_in_cell;
    
    % F-relaxation
    if strcmp( scheme, 'F' )

        % A CF interval is: C, F, ..., F, with at most m-1 F points. (All
        % intervals have m-1 F-points, except possibly the last, but this is 
        % handled by using the indptr array).
        for cpoint_idx = 1:nt_coarse % This loop is parallelizable.

            % We step from at the C point at the start of the current CF
            % interval, so get global index of C point (global among all points 
            % on the grid)
            start = cpoint_global_inds(cpoint_idx);     % Index of C-point at the start of the interval.
            stop  = cpoint_global_inds(cpoint_idx+1)-1; % Index of last F-point in the interval.

            % Step from t0 to t1
            for t0_idx = start:stop-1    
                t1_idx         = t0_idx + 1;
                
                step_status.t0_idx = t0_idx;
                step_status.t0     = t(t0_idx);
                step_status.dt     = t(t1_idx) - t(t0_idx);
                
                [utemp, MGRIT_object] = step(u_out_cell{t0_idx}, step_status, MGRIT_object);
                u_out_cell{t1_idx}    = utemp + g_cell{t1_idx};
            end
        end
    end

    % C-relaxation
    if strcmp( scheme, 'C' )

        % First C-point is only coupled to iself in the initial-value problem.
        u_out_cell{1} = g_cell{1};

        % Local index of C point (local among all C points)
        % This loop is parallelizable. Impoartantly, this is where we need 
        % to distinguish between the in and out array. Because the loop is 
        % implemented in serial, if the output array overwrites the input 
        % array and the coarsening factor is one, then C-relaxation will 
        % just do sequential time-stepping.
        for cpoint_idx = 2:nt_coarse 

            % Global index of C point (global among all points on the grid)
            cpoint_idx_global = cpoint_global_inds(cpoint_idx);

            % We step from the F point to the left of the current C-point
            t1_idx         = cpoint_idx_global;
            t0_idx         = t1_idx - 1;
            
            step_status.t0_idx = t0_idx;
            step_status.t0     = t(t0_idx);
            step_status.dt     = t(t1_idx) - t(t0_idx);
            
            [utemp, MGRIT_object] = step(u_in_cell{t0_idx}, step_status, MGRIT_object);
            u_out_cell{t1_idx}    = utemp + g_cell{t1_idx};
        end
    end        
end
% End of relaxation



%%
% Solve A(u) = b.
function [u_cell, MGRIT_object] = time_stepping_solve(b_cell, level, step, MGRIT_object)

    step_status.level = level;
    step_status.parent_function = 'time-stepping-solve';

    t      = MGRIT_object.hierarchy(level).t;
    u_cell = cell(size(b_cell));

    % Apply initial condition.
    u_cell{1} = b_cell{1};

    % Step over all other nt-1 points, from t(n) to t(n+1).
    for n = 1:numel(t)-1
        t0_idx = n;
        t1_idx = t0_idx + 1;
        
        step_status.t0_idx = n;
        step_status.t0     = t(t0_idx);
        step_status.dt     = t(t1_idx) - t(t0_idx);
        
        [utemp, MGRIT_object] = step(u_cell{t0_idx}, step_status, MGRIT_object);
        u_cell{t1_idx}        = utemp + b_cell{t1_idx};
    end
end
% End of time_stepping_solve


%% Residual at C-points only
function [rc_cell, MGRIT_object] = Cpoint_residual(u_cell, g_cell, level, step, MGRIT_object)

    step_status.level  = level;
    step_status.parent_function = 'C-point-residual';

    t         = MGRIT_object.hierarchy(level).t;
    %nt_coarse = MGRIT_object.hierarchy(level+1).nt;
    nt_coarse = numel(MGRIT_object.hierarchy(level).cpoint_indptr)-1;
    cf        = MGRIT_object.hierarchy(level).cf;

    % Number of C-points on current level.
    rc_cell = cell(nt_coarse, 1);

    % First C-point just couples to the initial conditon.
    rc_idx = 1;
    rc_cell{rc_idx} = g_cell{1} - u_cell{1};

    % Every other C-point is coupled to the F-point to its left.
    for CF_interval_idx = 2:nt_coarse
        rc_idx = rc_idx+1;

        % Last F point in the CF interval to left of C-point in question
        t0_idx = (CF_interval_idx-1) * cf;
        t1_idx = t0_idx + 1;

        step_status.t0_idx = t0_idx;
        step_status.t0     = t(t0_idx);
        step_status.dt     = t(t1_idx) - t(t0_idx);

        [uctemp, MGRIT_object] = step(u_cell{t0_idx}, step_status, MGRIT_object);
        rc_cell{rc_idx} = g_cell{t1_idx} + uctemp - u_cell{t1_idx};
    end

end
% End of Cpoint_residual

% % Compute r = b - A(u).
% function [r_cell, MGRIT_object] = full_residual(u_cell, b_cell, level, step, MGRIT_object)
% 
%     step_status.level = level;
%     step_status.parent_function = 'full-residual';
% 
%     t      = MGRIT_object.hierarchy(level).t;
%     r_cell = cell(size(b_cell));
% 
%     % Apply initial condition.
%     r_cell{1} = b_cell{1} - u_cell{1};
% 
%     % Step over all other nt-1 points, from t(n) to t(n+1).
%     for n = 1:numel(t)-1
%         t0_idx = n;
%         t1_idx = t0_idx + 1;
%         
%         step_status.t0_idx = n;
%         step_status.t0     = t(t0_idx);
%         step_status.dt     = t(t1_idx) - t(t0_idx);
%         
%         [utemp, MGRIT_object] = step(u_cell{t0_idx}, step_status, MGRIT_object);
%         r_cell{t1_idx}        = b_cell{t1_idx} + utemp - u_cell{t1_idx};
%     end
% end
% % End of full_residual