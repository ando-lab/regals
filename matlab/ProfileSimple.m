classdef ProfileSimple < ProfileClass
   properties
        % DEFINED BY SUPERCLASS
        % q

    end
    properties(Dependent = true)

        % INTERFACE REQUIRED BY SUPERCLASS (Nw,F,w,k,l,A)
        Nw % = Nq
        F  % sparse identity matrix
        w  % = q
        k  % = Nq
        L  % sparse identity matrix
        A  % error-weighted sparse identity matrix
        u0 % initial value

        % DEFINED BY SUPERCLASS
        % Nq  % = length(q)
        
        maxinfo
    end

    methods
        function obj = ProfileSimple(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end

        % GET METHODS (adapted from Humpty.m)
        function val = get.F(obj)
            val = speye(obj.Nq,obj.Nq);
        end
        function val = get.k(obj)
            val = obj.Nq;
        end
        function val = get.u0(obj)
            val = ones(obj.k,1);
        end
        function val = get.w(obj)
            val = obj.q(:);
        end
        function val = get.Nw(obj)
            val = obj.Nq;
        end
        function val = get.L(obj)
            val = obj.F;
        end
        function val = get.A(obj)
            val = obj.F;
            
        end
        function val = norm(~,u)
            val = sqrt(mean(u(:).*u(:)));
        end
        
        function val = get.maxinfo(obj)
            % estimate the maximum number of good parameters that can be
            % extracted from the data using this parameterization.

            val = numel(obj.q); % number of data points
        end

    end
end
