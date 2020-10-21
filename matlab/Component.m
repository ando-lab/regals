classdef Component
    properties
        Peak
        Profile
    end
    methods
        function obj = Component(Peak,Profile)
            obj.Peak = Peak;
            obj.Profile = Profile;
        end
    end

end
