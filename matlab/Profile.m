classdef Profile
    %PROFILE describes a 1D function in the q direction
    properties(SetAccess = immutable, GetAccess = private)
        regularizer
    end
    
    properties(SetAccess = private)
        type
        A    % cache value of obj.regularizer.A;
        L    % cache value of obj.regularizer.L;
        k    % cache value of obj.regularizer.k;
        u0   % cache value of obj.regularizer.u0;
        y0   % cache value of obj.regularizer.y0;
        Nq   % cache value of obj.regularizer.Nq;
        w    % cache value of obj.regularizer.w;
        maxinfo % cache value of obj.regularizer.maxinfo;
    end
    
    methods
        function obj = Profile(q,type,varargin)
            obj.type = lower(type);
            switch obj.type
                case {'','simple'}
                    obj.regularizer = ProfileSimple('q',q,varargin{:});
                case 'smooth'
                    obj.regularizer = ProfileSmooth('q',q,varargin{:});
                case 'realspace'
                    obj.regularizer = ProfileRealSpace('q',q,varargin{:});
                otherwise
                    error('unknown profile regularizer');
            end
            
            obj.A = obj.regularizer.A;
            obj.L = obj.regularizer.L;
            obj.k = obj.regularizer.k;
            obj.u0 = obj.regularizer.u0;
            obj.y0 = obj.regularizer.y0;
            obj.Nq = obj.regularizer.Nq;
            obj.w = obj.regularizer.w;
            obj.maxinfo = obj.regularizer.maxinfo;
        end
        
        function val = norm(obj,u)
            val = obj.regularizer.norm(u);
        end
        
    end
end

