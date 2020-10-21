classdef ProfileRealSpace < ProfileClass
    properties
        % INTERFACE REQUIRED BY SUPERCLASS (Nw)
        Nw = 50 % default number of samples

        % DEFINED BY SUPERCLASS
        % q

        % REQUIRED PROPERTIES
        dmax
        isZeroAtDmax = true
        isZeroAtR0 = true

    end

    properties(Dependent = true)

        % INTERFACE REQUIRED BY SUPERCLASS (F,w,k,l,A)
        F  % fourier transform operator: I(q) = F*P(r) [here, P(r) is of length Nr]
        w  % r vector
        k  % number of coefficients of P(r) which are free
        L  % second derivative operator: P''(r) = Lr*P(r) [here, P(r) is of length kr]
        A  % truncated and error-weighted Fourier transform operator (Nq by kr)
        u0 % an initial value for parameters

        % DEFINED BY SUPERCLASS:
        % Nq % = length(q)
    end
    properties(Dependent = true, Access = private)
        dw % spacing between samples in w vector
    end
    methods
        function obj = ProfileRealSpace(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end

        % SAXS PROFILE-RELATED GET METHODS (see Nifty.m)
        function val = get.F(obj)
            val = 4*pi*obj.dw*sinc(obj.q(:)*obj.w'/pi);
            val(:,[1,end]) = 0.5*val(:,[1,end]); % trapezoid rule
        end
        function val = get.k(obj)
            val = obj.Nw - obj.isZeroAtDmax - obj.isZeroAtR0;
        end
        function val = get.dw(obj)
            val = obj.dmax/(obj.Nw-1);
        end
        function val = get.w(obj)
            val = obj.dmax*linspace(0,1,obj.Nw)';
        end
        function val = get.u0(obj)
            val = 1 - ( 2*obj.w/obj.dmax - 1).^2;
            val = val/(4*pi*obj.dw*sum(val)); % normalize to I(0)=1
            if obj.isZeroAtDmax
                val = val(1:(end-1));
            end
            if obj.isZeroAtR0
                val = val(2:end);
            end
        end
        function val = norm(obj,u)
            % return I(0)
            weight = 4*pi*obj.dw*ones(obj.Nw,1);
            weight([1,end]) = 0.5*weight([1,end]); % trapezoid rule
            if obj.isZeroAtDmax
                weight = weight(1:(end-1));
            end
            if obj.isZeroAtR0
                weight = weight(2:end);
            end
            val = sum(weight.*u(:));
        end
        %function val = get.Nq(obj)
        %    val = length(obj.q);
        %end
        function val = get.L(obj)
            val = sparse(1:(obj.Nw-2),1:(obj.Nw-2),-.5,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),2:(obj.Nw-1),1,obj.Nw-2,obj.Nw) + ...
                sparse(1:(obj.Nw-2),3:(obj.Nw),-.5,obj.Nw-2,obj.Nw);
            if obj.isZeroAtDmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtR0
                val = val(:,2:end);
            end
        end
        function val = get.A(obj)
            val = obj.F;

            if obj.isZeroAtDmax
                val = val(:,1:(end-1));
            end
            if obj.isZeroAtR0
                val = val(:,2:end);
            end
        end
    end
end
