classdef ProfileSmooth < ProfileClass
   properties
        % INTERFACE REQUIRED BY SUPERCLASS (Nw)
        Nw = 50 % default number of samples

        % DEFINED BY SUPERCLASS
        % q
    end
    properties(Dependent = true)

        % INTERFACE REQUIRED BY SUPERCLASS (F,w,k,l,A)
        F  % linear interpolation operator: I(q) = F*ws
        w  % vector of points sampling the peak
        k  % number of coefficients of ws which are free
        L  % second derivative operator: W''(s) = Ls*W(s) [here, W(s) is of length ks]
        A  % truncated and error-weighted interpolation matrix (Ns by ks)
        u0 % inital constant value

        % DEFINED BY SUPERCLASS
        % Nq  % = length(q)
    end
    properties(Dependent = true, Access = private)
        dw % spacing between samples in w vector
    end
    properties(Dependent = true, Access = private)
        qmin
        qmax
    end

    methods
        function obj = ProfileSmooth(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end

        % GET METHODS
        function val = get.F(obj)
            thisq = obj.q(:);
            [~,ix] = histc(thisq,[-Inf;obj.w;Inf]);
            ix = ix-1;
            ix(ix==0) = 1;
            ix(ix==obj.Nw) = obj.Nw-1;
            v = (1:obj.Nq)';
            u = (thisq - obj.w(ix))*(1/obj.dw);
            val = sparse([v;v],[ix;ix+1],[1-u;u],obj.Nq,obj.Nw);
        end
        function val = get.k(obj)
            val = obj.Nw;
        end
        function val = get.u0(obj)
            val = ones(obj.k,1);
        end
        function val = get.qmin(obj)
            val = min(obj.q);
        end
        function val = get.qmax(obj)
            val = max(obj.q);
        end

        function val = get.dw(obj)
            val = (obj.qmax - obj.qmin)/(obj.Nw-1);
        end
        function val = get.w(obj)
            val = linspace(obj.qmin,obj.qmax,obj.Nw)';
        end
        %function val = get.Nq(obj) % <--- defined in superclass
        %    val = length(obj.q);
        %end
        function val = get.L(obj)
            val = sparse(1:(obj.Nw-2),1:(obj.Nw-2),-.5,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),2:(obj.Nw-1),1,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),3:(obj.Nw),-.5,obj.Nw-2,obj.Nw);
        end
        function val = get.A(obj)
            val = obj.F;

        end
        function val = norm(~,u)
            val = sqrt(mean(u(:).*u(:)));
        end

    end
end
