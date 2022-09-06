classdef Nifty
    properties
        q
        int
        err
        dmax
        Nr = 50
        isZeroAtDmax = true
        isZeroAtZero = true
    end
    properties(Dependent = true)
        F
        r
        dr
        k
        L
        A
        Nq
        b
    end
    methods
        function obj = Nifty()
        end
        function val = get.F(obj)
            val = 4*pi*obj.dr*sinc(obj.q(:)*obj.r'/pi);
            val(:,[1,end]) = 0.5*val(:,[1,end]); % trapezoid rule
        end
        function val = get.k(obj)
            val = obj.Nr - obj.isZeroAtDmax - obj.isZeroAtZero;
        end
        function val = get.dr(obj)
            val = obj.dmax/(obj.Nr-1);
        end
        function val = get.r(obj)
            val = obj.dmax*linspace(0,1,obj.Nr)';
        end
        function val = get.Nq(obj)
            val = length(obj.q);
        end
        function val = get.L(obj)
            val = sparse(1:(obj.Nr-2),1:(obj.Nr-2),-.5,obj.Nr-2,obj.Nr) + ...
                sparse(1:(obj.Nr-2),2:(obj.Nr-1),1,obj.Nr-2,obj.Nr) + ...
                sparse(1:(obj.Nr-2),3:(obj.Nr),-.5,obj.Nr-2,obj.Nr);
            if obj.isZeroAtDmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtZero
                val = val(:,2:end);
            end
        end
        function val = get.A(obj)
            val = obj.F;
            if length(obj.err)==1
                val = val/obj.err;
            elseif isempty(obj.err)
                % do nothing
            else
                val = repmat(1./obj.err(:),1,obj.Nr).*val;
            end
            if obj.isZeroAtDmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtZero
                val = val(:,2:end);
            end
        end
        function val = get.b(obj)
            val = obj.int(:)./obj.err(:);
        end
        function [w,alpha] = ift(obj,alpha)
            Lmat = obj.L;
            Amat = obj.A;
            AA = Amat'*Amat;
            LL = full(Lmat'*Lmat);
            if isempty(alpha)
                alpha = trace(AA)/trace(LL);
            end
            w = (AA + alpha*LL)\(Amat'*obj.b);
            if obj.isZeroAtDmax
                w = [w;0];
            end
            if obj.isZeroAtZero
                w = [0;w];
            end
        end

        function [y,y0,w0] = make_test_data(obj,sigma)
            if nargin==1 || isempty(sigma)
                sigma = obj.err;
            end
            x = obj.r/obj.dmax;
            w0 = 1-(2*x-1).^2;
            w0 = w0.*x.*exp(-2*x);
            w0 = w0 / (4*pi*obj.dr*sum(w0)); % normalize (I(0) = 1)
            y0 = obj.F*w0;
            if ~isempty(sigma)
                y = y0  + randn([obj.Nq,1]).*sigma;
            else
                y = [];
            end
        end
    end
        
    methods(Static)
        function [r,w,ireg,alpha] = sprite(q,int,err,varargin)
                         
            N = Nifty();
            N.q = q;
            N.int = int;
            N.err = err;
            
            % assign optional arguments
            for j=1:2:length(varargin)
                N.(varargin{j}) = varargin{j+1};
            end
            
            [w,alpha] = N.ift([]);
            r = N.r;
            ireg = N.F*w;
        end
    end
end

