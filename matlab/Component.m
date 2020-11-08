classdef Component
    properties
        Concentration
        Profile
    end
    methods
        function obj = Component(Concentration,Profile)
            obj.Concentration = Concentration;
            obj.Profile = Profile;
        end
    end

end
