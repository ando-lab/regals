classdef ConcentrationSmooth < ConcentrationClass
    properties
        
        % INTERFACE REQUIRED BY SUPERCLASS (Nw)
        Nw = 50 % default number of samples
        
        % DEFINED BY SUPERCLASS
        % x
        
        % REQUIRED PROPERTIES
        xmin
        xmax
        isZeroAtXmin = true
        isZeroAtXmax = true
        
    end
    
    properties(Dependent = true)
        % INTERFACE REQUIRED BY SUPERCLASS (F,w,k,l,A)
        F  % linear interpolation operator: C(x) = I*wx
        w  % vector of points sampling the concentration
        k  % number of coefficients of w which are free
        L  % second derivative operator: W''(x) = Lx*W(x) [here, W(x) is of length k]
        A  % truncated and error-weighted interpolation matrix (Nx by k)
        u0 % initial value of u to choose
        % DEFINED BY SUPERCLASS:
        % Nx % = length(x)
        maxinfo
    end
    properties(Dependent = true, Access = private)
        dw % spacing between samples in w vector
    end
    methods
        function obj = ConcentrationSmooth(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end
        
        % CONCENTRATION PROFILE-RELATED GET METHODS (see Humpty.m)
        function val = get.F(obj)
            %             thisx = obj.x(:);
            %             [~,ix] = histc(thisx,[-Inf;obj.w]);
            %             ix = ix-1;
            %             ix(ix==obj.Nw) = obj.Nw - 1; % if equal to last bin edge, assign to last bin
            %             isInConcentration = ix > 0 & ix < obj.Nw;
            %             v = (1:obj.Nx)';
            %             ix = ix(isInConcentration);
            %             v = v(isInConcentration);
            %            u = (thisx(isInConcentration) - obj.w(ix))*(1/obj.dw);
            [ix,v] = obj.xindex();
            thisx = obj.x(:);
            u = (thisx(v) - obj.w(ix))*(1/obj.dw);
            val = sparse([v;v],[ix;ix+1],[1-u;u],obj.Nx,obj.Nw);
        end
        
        function val = get.maxinfo(obj)
            % estimate the maximum number of good parameters that can be
            % extracted from the data using this parameterization.
            ix = obj.xindex();
            Nd = numel(ix); % number of data points in range
            
            val = min([Nd,obj.k]); % Whichever is less: number of free control points or number of data points between xmin, xmax
            
            % note, this does not consider what happens when a data point
            % is on the boundary, and the boundary condition is set to zero
            % (in that case, the data point no longer contributes)
        end
        
        function val = get.k(obj)
            val = obj.Nw - obj.isZeroAtXmin - obj.isZeroAtXmax;
        end
        function val = get.dw(obj)
            val = (obj.xmax - obj.xmin)/(obj.Nw-1);
        end
        function val = get.w(obj)
            val = linspace(obj.xmin,obj.xmax,obj.Nw)';
        end
        function val = get.u0(obj)
            if obj.isZeroAtXmax && obj.isZeroAtXmin
                val = 1 - ( 2*(obj.w-obj.xmin)/(obj.xmax-obj.xmin) - 1).^2;
            elseif obj.isZeroAtXmax && ~obj.isZeroAtXmin
                val = (obj.xmax-obj.w)/(obj.xmax-obj.xmin);
            elseif ~obj.isZeroAtXmax && obj.isZeroAtXmin
                val = (obj.w-obj.xmin)/(obj.xmax-obj.xmin);
            else
                val = ones(obj.Nw,1);
            end
            
            if obj.isZeroAtXmax
                val = val(1:(end-1));
            end
            if obj.isZeroAtXmin
                val = val(2:end);
            end
        end
        function val = get.L(obj)
            val = sparse(1:(obj.Nw-2),1:(obj.Nw-2),-.5,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),2:(obj.Nw-1),1,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),3:(obj.Nw),-.5,obj.Nw-2,obj.Nw);
            if obj.isZeroAtXmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtXmin
                val = val(:,2:end);
            end
        end
        function val = get.A(obj)
            val = obj.F;
            
            if obj.isZeroAtXmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtXmin
                val = val(:,2:end);
            end
        end
    end
    methods(Access = private)
        function [ix,v] = xindex(obj)
            % helper function for get.F
            thisx = obj.x(:);
            [~,ix] = histc(thisx,[-Inf;obj.w]);
            ix = ix-1;
            ix(ix==obj.Nw) = obj.Nw - 1; % if equal to last bin edge, assign to last bin
            isInConcentration = ix > 0 & ix < obj.Nw;
            v = (1:obj.Nx)';
            ix = ix(isInConcentration);
            v = v(isInConcentration);
        end
    end
end
