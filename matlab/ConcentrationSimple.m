classdef ConcentrationSimple < ConcentrationClass
    properties
        % DEFINED BY SUPERCLASS
        % x

        % REQUIRED PROPERTIES
        xmin
        xmax
    end

    properties(Dependent = true)
        % INTERFACE REQUIRED BY SUPERCLASS (F,w,k,l,A)
        Nw % number of points in x between xmin and xmax inclusive
        F  % sparse identity matrix (Nx by Nw)
        w  % points in x between xmin and xmax inclusive
        k  % = Nw
        L  % sparse identity matrix (Nw by Nw)
        A  % error-weighted identity matrix (Nx by Nw)
        u0 % initial value of w to choose

        % DEFINED BY SUPERCLASS:
        % Nx % = length(x)
    end
    methods
        function obj = ConcentrationSimple(varargin)
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end

        % CONCENTRATION PROFILE-RELATED GET METHODS (see Humpty.m)
        function val = get.F(obj)
            thisx = obj.x(:);
            thisw = obj.w;
            [isInConcentration,indInW] = ismember(thisx,thisw);
            v = (1:obj.Nx)';
            val = sparse(v(isInConcentration),indInW(isInConcentration),1,obj.Nx,length(thisw));
        end
        function val = get.k(obj)
            val = obj.Nw;
        end
        function val = get.w(obj)
            thisx = obj.x(:);
            val = unique(thisx(thisx >= obj.xmin & thisx <= obj.xmax));
        end
        function val = get.Nw(obj)
            val = length(obj.w);
        end
        function val = get.u0(obj)
            val = ones(obj.k,1);
        end
        function val = get.L(obj)
            val = speye(obj.Nw);
        end
        function val = get.A(obj)
            val = obj.F;
        end

    end
end
