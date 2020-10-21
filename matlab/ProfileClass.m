classdef (Abstract) ProfileClass
    properties
        q
    end
    properties (Dependent = true)
        Nq
        y0
    end
    properties (Abstract)
        Nw
        w
        k
        F
        L
        A
        u0
    end
    methods
        function val = get.Nq(obj)
            val = length(obj.q);
        end
        function val = get.y0(obj)
            val = obj.A*obj.u0;
        end
    end
end