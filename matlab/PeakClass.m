classdef (Abstract) PeakClass
    properties
        x
    end
    properties (Dependent = true)
        Nx
        y0
    end
    properties (Abstract)
        Nw
        w
        k
        F
        L
        A
        u0 % initial value of w to choose
    end
    methods
        function val = get.Nx(obj)
            val = length(obj.x);
        end
        function val = get.y0(obj)
            val = obj.A*obj.u0;
        end
        function val = norm(~,u)
            val = sqrt(mean(u(:).*u(:)));
        end
    end
end