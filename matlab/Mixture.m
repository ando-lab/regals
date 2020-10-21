classdef Mixture
    % MIXTURE - class representing the mixture model and functions for
    % fitting either peaks or profiles.
    %
    % This is used by algorithms such as REGALS to store the current state
    % of the optimization, and to prepare the matrices for least-squares.
    properties
        Components
        peakLambda
        profileLambda
        uPeak
        uProfile
    end

    properties(Dependent = true)
        Nc % number of components
        Nx
        Nq
        kPeak
        kProfile
        peaks
        profiles
        peakNorm
        profileNorm
        peakH
        profileH
        Ireg
    end

    methods
        function obj = Mixture(Comp,varargin)

            obj.Components = Comp;

            % assign arguments
            for j=1:2:(numel(varargin)-1)
                obj.(varargin{j}) = varargin{j+1};
            end

            if isempty(obj.uPeak)
                % get default values for uPeak
                obj.uPeak = cell(obj.Nc,1);
                for j=1:obj.Nc
                    obj.uPeak{j} = obj.Components(j).Peak.u0(:);
                end
            end

            if isempty(obj.uProfile)
                % get default values for uProfile
                obj.uProfile = cell(obj.Nc,1);
                for j=1:obj.Nc
                    obj.uProfile{j} = obj.Components(j).Profile.u0(:);
                end
            end

            if isempty(obj.peakLambda)
                obj.peakLambda = zeros(1,obj.Nc);
            end

            if isempty(obj.profileLambda)
                obj.profileLambda = zeros(1,obj.Nc);
            end
        end

        function val = get.Nc(obj)
            val = numel(obj.Components);
        end

        function val = get.Nx(obj)
            val = obj.Components(1).Peak.Nx;
        end

        function val = get.Nq(obj)
            val = obj.Components(1).Profile.Nq;
        end

        function k = get.kPeak(obj)
            Peaks = [obj.Components.Peak];
            k = [Peaks.k];
        end

        function k = get.kProfile(obj)
            Profiles = [obj.Components.Profile];
            k = [Profiles.k];
        end

        function I = get.Ireg(obj)
            I = obj.profiles*obj.peaks';
        end

        function ll = estimatePeakLambda(obj,err,ng)
            AA = obj.peakProblem(sparse(obj.Nq,obj.Nx),err);

            % break into sub-matrices for each component
            AA = mat2cell(AA,obj.kPeak,obj.kPeak);

            ll = zeros(1,obj.Nc);
            for k=1:obj.Nc
                Lk = obj.Components(k).Peak.L;
                [~,d] = eig(full(AA{k,k}),full(Lk'*Lk));
                ll(k) = ng2lambda(diag(d),ng(k));
            end
        end

        function ll = estimateProfileLambda(obj,err,ng)
            AA = obj.profileProblem(sparse(obj.Nq,obj.Nx),err);

            % break into sub-matrices for each component
            AA = mat2cell(AA,obj.kProfile,obj.kProfile);

            ll = zeros(1,obj.Nc);
            for k=1:obj.Nc
                Lk = obj.Components(k).Profile.L;
                [~,d] = eig(full(AA{k,k}),full(Lk'*Lk));
                ll(k) = ng2lambda(diag(d),ng(k));
            end
        end

        function [AA,Ab] = profileProblem(obj,I,err)
            % assumes equal errorbars for each q-bin

            w = 1./mean(err,2); % average over columns

            nc = numel(obj.Components);
            AA = cell(nc,nc);
            Ab = cell(nc,1);

            Profiles = [obj.Components.Profile];
            A = {Profiles.A};

            % apply error weight to A-matrices
            for k = 1:nc
                A{k} = w.*A{k};
            end

            % apply error weight to data matrix
            D = w.*I;

            c = obj.peaks;

            % calculate AA = A'A
            for k1 = 1:nc
                for k2 = 1:nc
                    AA{k1,k2} = (c(:,k1)'*c(:,k2))*(A{k1}'*A{k2});
                end
            end

            % calculate rhs: Ab = A'b
            for k = 1:nc
                Ab{k} = A{k}'*(D*c(:,k));
            end

            % convert from cell to array
            AA = cell2mat(AA);
            Ab = cell2mat(Ab);
        end

        function [AA,Ab] = peakProblem(obj,I,err)
            % assumes equal errorbars for each q-bin

            w = 1./mean(err,2); % average over q-bins

            nc = numel(obj.Components);
            AA = cell(nc,nc);
            Ab = cell(nc,1);

            Peaks = [obj.Components.Peak];
            A = {Peaks.A};

            % apply error weight to data matrix
            D = (w.*I); % fit the transpose

            y = obj.profiles;

            % apply error weight to each profile
            y = w.*y;

            % calculate AA = A'A
            for k1 = 1:nc
                for k2 = 1:nc
                    AA{k1,k2} = (y(:,k1)'*y(:,k2))*(A{k1}'*A{k2});
                end
            end

            % calculate rhs: Ab = A'b
            for k = 1:nc
                Ab{k} = A{k}'*(D'*y(:,k));
            end

            % convert from cell to array
            AA = cell2mat(AA);
            Ab = cell2mat(Ab);
        end

        function [Ik,sigmak] = extractProfile(obj,I,err,k)
            notk = setdiff(1:obj.Nc,k);
            c = obj.peaks;
            y = obj.profiles;
            D = I - y(:,notk)*c(:,notk)';
            ck = c(:,k);
            m = ck/(ck'*ck); % coefficients in the sum
            Ik = D*m;
            sigmak = sqrt(err.^2*m.^2);
        end

        function [pk,sigmak] = extractPeak(obj,I,err,k)
            notk = setdiff(1:obj.Nc,k);
            c = obj.peaks;
            y = obj.profiles;
            D = I - y(:,notk)*c(:,notk)';
            yk = y(:,k);
            w = 1./mean(err,2); % error weights
            m = (w.^2.*yk)/(yk'*(w.^2.*yk)); % coefficients in the sum
            pk = D'*m;
            sigmak = sqrt(err'.^2*m.^2);
        end

        function H = get.peakH(obj)
            B = cell(1,obj.Nc);
            for j=1:obj.Nc
                L = obj.Components(j).Peak.L;
                B{j} = L*sqrt(obj.peakLambda(j));
            end
            B = blkdiag(B{:});
            H = B'*B;
        end

        function H = get.profileH(obj)
            B = cell(1,obj.Nc);
            for j=1:obj.Nc
                L = obj.Components(j).Profile.L;
                B{j} = L*sqrt(obj.profileLambda(j));
            end
            B = blkdiag(B{:});
            H = B'*B;
        end

        function p = get.profiles(obj)
            Profiles = [obj.Components.Profile];
            p = zeros(obj.Nq,obj.Nc);
            for j=1:obj.Nc
                p(:,j) = Profiles(j).A*obj.uProfile{j};
            end
        end

        function p = get.peaks(obj)
            Peaks = [obj.Components.Peak];
            p = zeros(obj.Nx,obj.Nc);
            for j=1:obj.Nc
                p(:,j) = Peaks(j).A*obj.uPeak{j};
            end
        end

        function n = get.peakNorm(obj)
            n = zeros(1,obj.Nc);
            for j=1:obj.Nc
                n(j) = obj.Components(j).Peak.norm(obj.uPeak{j});
            end
        end

        function n = get.profileNorm(obj)
            n = zeros(1,obj.Nc);
            for j=1:obj.Nc
                n(j) = obj.Components(j).Profile.norm(obj.uProfile{j});
            end
        end

    end
end


function lambda = ng2lambda(dd,ng)

ng0 = nnz(isinf(dd)); % minimum number of free parameters

dd = dd(dd>=0 & ~isinf(dd)); % remove negative and infinite values

lambdaList = logspace(max(log10(dd)) + 2,min(log10(dd)) - 2,51);

% invert by linear interpolation on a log scale
ngList = zeros(size(lambdaList));
for j=1:numel(lambdaList)
    ngList(j) = ng0 + sum(dd./(dd + lambdaList(j)));
end

if ng < ngList(1) % edge case 1
    lambda = Inf;
elseif ng > ngList(end) % edge case 2
    lambda = 0;
else % interpolation
    lambda = 10^interp1q(ngList',log10(lambdaList)',ng);
end

end
