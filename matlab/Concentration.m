classdef Concentration
    %CONCENTRATION describes a 1D function in the x direction
    properties(SetAccess = immutable, GetAccess = private)
        regularizer
    end
    
    properties(SetAccess = private)
        type % simple or smooth
        A    % cache obj.regularizer.A;
        L    % cache obj.regularizer.L;
        k    % cache obj.regularizer.k;
        u0   % cache obj.regularizer.u0;
        y0   % cache obj.regularizer.y0;
        Nx   % cache obj.regularizer.Nx;
        w    % cache obj.regularizer.w;
        maxinfo % cache obj.regularizer.maxinfo
    end
    
    methods
        function obj = Concentration(x,type,varargin)
            obj.type = lower(type);
            switch obj.type
                case {'','simple'}
                    obj.regularizer = ConcentrationSimple('x',x,varargin{:});
                case 'smooth'
                    obj.regularizer = ConcentrationSmooth('x',x,varargin{:});
                otherwise
                    error('unknown concentration regularizer');
            end
            
            obj.A = obj.regularizer.A;
            obj.L = obj.regularizer.L;
            obj.k = obj.regularizer.k;
            obj.u0 = obj.regularizer.u0;
            obj.y0 = obj.regularizer.y0;
            obj.Nx = obj.regularizer.Nx;
            obj.w = obj.regularizer.w;
            obj.maxinfo = obj.regularizer.maxinfo;
            
        end
        
        function val = norm(obj,u)
            val = obj.regularizer.norm(u);
        end
        
    end
    
end

