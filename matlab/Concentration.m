classdef Concentration
    %CONCENTRATION describes a 1D function in the x direction
    properties(SetAccess = immutable, GetAccess = private)
        regularizer
    end
    
    properties(SetAccess = private)
        A    % cache obj.regularizer.A;
        L    % cache obj.regularizer.L;
        k    % cache obj.regularizer.k;
        u0   % cache obj.regularizer.u0;
        y0   % cache obj.regularizer.y0;
        Nx   % cache obj.regularizer.Nx;
        w    % cache obj.regularizer.w;
    end
    
    methods
        function obj = Concentration(x,type,varargin)
            if nargin > 1
                switch lower(type)
                    case {'','simple'}
                        obj.regularizer = ConcentrationSimple('x',x,varargin{:});
                    case 'smooth'
                        obj.regularizer = ConcentrationSmooth('x',x,varargin{:});
                    otherwise
                        error('unknown concentration regularizer');
                end
            end
            
            obj.A = obj.regularizer.A;
            obj.L = obj.regularizer.L;
            obj.k = obj.regularizer.k;
            obj.u0 = obj.regularizer.u0;
            obj.y0 = obj.regularizer.y0;
            obj.Nx = obj.regularizer.Nx;
            obj.w = obj.regularizer.w;
            
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
        function val = get.Nx(obj)
            val = obj.regularizer.Nx;
        end
        function val = norm(obj,u)
            val = obj.regularizer.norm(u);
        end
    end
    
end

