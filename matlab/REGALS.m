classdef REGALS
    %REGALS - regularized alternating least squares
    %
    % May 15, 2020 - commented out experimental functions
    properties
        I
        err
    end

    methods
        function obj = REGALS(I,err,varargin)
            %REGALS
            obj.I = I;
            obj.err = err;
            if ~isempty(varargin)
                for j=1:2:length(varargin)
                    obj.(varargin{j}) = varargin{j+1};
                end
            end
        end

        function Mix = fitPeaks(obj,Mix)
            %FITPEAKS

            % calculate matrices
            H = Mix.peakH;
            [AA,Ab] = Mix.peakProblem(obj.I,obj.err);

            % least-squares fit
            u = (AA + H)\Ab;

            % repartition u
            u = mat2cell(u,Mix.kPeak,1);

            % assign to output
            Mix.uPeak = u;
        end

        function Mix = fitProfiles(obj,Mix)
            %FITPROFILES

            % calculate matrices
            H = Mix.profileH;
            [AA,Ab] = Mix.profileProblem(obj.I,obj.err);

            % least-squares fit
            u = (AA + H)\Ab;

            % repartition u
            u = mat2cell(u,Mix.kProfile,1);

            % normalize
            n = Mix.profileNorm;
            for j=1:Mix.Nc
                u{j} = u{j}/n(j);
            end

            % assign to output
            Mix.uProfile = u;
        end

        function [NewMix,params,resid] = step(obj,Mix)
            % step - take a REGALS step

            % fit profiles then peaks
            NewMix = obj.fitPeaks(obj.fitProfiles(Mix));

            % calculate goodness of fit
            resid = (obj.I - NewMix.Ireg)./obj.err;
            x2 = mean(resid(:).^2);

            % calculate deltas
            delta_peak = sum(abs(NewMix.peaks - Mix.peaks),1);
            delta_profile = sum(abs(NewMix.profiles - Mix.profiles),1);
            delta_uPeak = cellfun(@(v1,v2) sum(abs(v1-v2)),...
                            NewMix.uPeak,Mix.uPeak);
            delta_uProfile = cellfun(@(v1,v2) sum(abs(v1-v2)),...
                            NewMix.uProfile,Mix.uProfile);

            % return various parameters
            params = struct(...
                'x2',x2,...
                'delta_peak',delta_peak,...
                'delta_profile',delta_profile,...
                'delta_uPeak',delta_uPeak',... % <- transpose to row vector
                'delta_uProfile',delta_uProfile');
        end

        function [Mix,params,resid,exitCond] = run(obj,Mix,stopFun,updateFun)
            % run - generalized function for iterative REGALS
            %
            % [Mix,params,resid,exitCond] = run(Mix,stopFun,updateFun)
            %
            % stopFun - function handle that returns true if stop conditions
            % are met. The function is called as:
            %
            %   [tf,exitCond] = stopFun(iter,params)
            %
            % the default is:
            %    stopFun = @(iter,params) deal(iter <  10,'maxIter')
            %
            % updateFun - optional function for displaying or logging output
            % while function is evaluating. The function is called as:
            %
            %   updateFun(iter,NewMix,params,resid)
            %
            % (for more information on arguments, see REGALS.step).

            % set default functions
            if nargin < 3 || isempty(stopFun)
                stopFun = @(iter,params) deal(iter >=  10,'maxIter');
            end
            if nargin < 4 || isempty(updateFun)
                updateFun = @(iter,NewMix,params,resid) true; % do nothing
            end

            iter = 0;
            while true
                iter = iter + 1;
                [Mix,params,resid] = obj.step(Mix);

                updateFun(iter,Mix,params,resid);

                % check whether or not to stop
                [tf,exitCond] = stopFun(iter,params);
                if tf
                    break;
                end
            end

        end

    end

end
