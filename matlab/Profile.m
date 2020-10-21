classdef Profile
    %PROFILE describes a 1D function in the q direction
    properties
        regularizer
    end
    
    properties(Dependent = true)
        A    % wrapper function for obj.regularizer.A;
        L    % wrapper function for obj.regularizer.L;
        k    % wrapper function for obj.regularizer.k;
        u0   % wrapper function for obj.regularizer.u0;
        y0   % wrapper function for obj.regularizer.y0;
        Nq   % wrapper function for obj.regularizer.Nq;
    end
    
    methods
        function obj = Profile(q,type,varargin)
            if nargin > 1
                switch lower(type)
                    case {'','simple'}
                        obj.regularizer = ProfileSimple('q',q,varargin{:});
                    case 'smooth'
                        obj.regularizer = ProfileSmooth('q',q,varargin{:});
                    case 'realspace'
                        obj.regularizer = ProfileRealSpace('q',q,varargin{:});
                    otherwise
                        error('unknown profile regularizer');
                end
            end
        end
        
        function val = get.A(obj)
            val = obj.regularizer.A;
        end
        function val = get.L(obj)
            val = obj.regularizer.L;
        end
        function val = get.k(obj)
            val = obj.regularizer.k;
        end
        function val = get.u0(obj)
            val = obj.regularizer.u0;
        end
        function val = get.y0(obj)
            val = obj.regularizer.y0;
        end
        function val = get.Nq(obj)
            val = obj.regularizer.Nq;
        end
        function val = norm(obj,u)
            val = obj.regularizer.norm(u);
        end
    end
end

