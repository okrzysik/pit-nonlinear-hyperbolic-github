% Class implementing the Burgers equation
classdef burgers < cons_law_scalar 

    properties

    end
    
    methods

        % Constructor
        function obj = burgers(pde_pa_)
            % Call parent's constructor
            obj = obj@cons_law_scalar(pde_pa_, true);
            
        end
        
        function f = flux(obj, u)
           f = 0.5 * u.^2;
        end
        
        function fprime = flux_jacobian(obj, u)
           fprime = u;
        end
        
        function id = id(obj)
            id = 'burgers';
        end
        
        function id = id_abbreviation(obj)
            id = 'B';
        end
        
        
        % Evaluate the initial condition. The below text is in reference to
        % the solution of Burgers equation.
        function u0 = initial_condition(obj, x)

            switch obj.pde_pa.ic_id
                % moving shock    
                case 1
                    u0 = 0.25 - 0.75*sin(pi*x);
                
                % Shock forms along the x = 0 line and dissipates away. See 
                % same kind of problem in Hesthaven (2017), Figs. 1.3 & 1.4.
                case 2    
                    u0 = -sin(pi*x); 

                    % Riemann problem. Right-moving shock and rarefaction wave    
                case 3
                    u0 = zeros(size(x));
                    I = find(x > -1/2 & x < 0);
                    u0(I) = 1;
                    
                case 4
                    u0 = 0.25 - 0.75*sin(pi*x) + 0.2*sin(12*pi*x);
                    
                case 5    
                    u0 = 0.5*(sin(2*pi*x) + sin(3*pi*x).^2);
                    
                case 6
                    u0 = 0.8 + 0.2*cos(pi*x);
                    
                otherwise
                    error('ic_id = %d not recognised', obj.ic_id)
                     
            end
        
        % End of initial_condition 
        end
        
    % End of methods.    
    end

% End of class.
end                    