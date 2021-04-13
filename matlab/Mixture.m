classdef Mixture
    % MIXTURE - class representing the mixture model and functions for
    % fitting either concentrations or profiles.
    %
    % This is used by algorithms such as REGALS to store the current state
    % of the optimization, and to prepare the matrices for least-squares.
    properties
        Components
        concentrationLambda
        profileLambda
        uConcentration
        uProfile
    end

    properties(Dependent = true)
        Nc % number of components
        Nx
        Nq
        kConcentration
        kProfile
        concentrations
        profiles
        concentrationNorm
        profileNorm
        concentrationH
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

            if isempty(obj.uConcentration)
                % get default values for uConcentration
                obj.uConcentration = cell(obj.Nc,1);
                for j=1:obj.Nc
                    obj.uConcentration{j} = obj.Components(j).Concentration.u0(:);
                end
            end

            if isempty(obj.uProfile)
                % get default values for uProfile
                obj.uProfile = cell(obj.Nc,1);
                for j=1:obj.Nc
                    obj.uProfile{j} = obj.Components(j).Profile.u0(:);
                end
            end

            if isempty(obj.concentrationLambda)
                obj.concentrationLambda = zeros(1,obj.Nc);
            end

            if isempty(obj.profileLambda)
                obj.profileLambda = zeros(1,obj.Nc);
            end
        end

        function val = get.Nc(obj)
            val = numel(obj.Components);
        end

        function val = get.Nx(obj)
            val = obj.Components(1).Concentration.Nx;
        end

        function val = get.Nq(obj)
            val = obj.Components(1).Profile.Nq;
        end

        function k = get.kConcentration(obj)
            Concentrations = [obj.Components.Concentration];
            k = [Concentrations.k];
        end

        function k = get.kProfile(obj)
            Profiles = [obj.Components.Profile];
            k = [Profiles.k];
        end

        function I = get.Ireg(obj)
            I = obj.profiles*obj.concentrations';
        end

        function ll = estimateConcentrationLambda(obj,err,ng)
            
            if nargin < 3 || isempty(ng) % need to estimate ng first
                ng = zeros(1,obj.Nc);
                for j=1:obj.Nc
                    C = obj.Components(j).Concentration;
                    switch C.type
                        case 'simple' % no regularization
                            ng(j) = Inf;
                        case 'smooth' % just a little smoothing
                            ng(j) = 0.8*C.maxinfo;
                        otherwise
                            error('concentration type not implemented?');
                    end
                end
                fprintf(1,'estimating concentration lambda with ng = %s\n',mat2str(ng));
            end
            
            AA = obj.concentrationProblem(sparse(obj.Nq,obj.Nx),err);

            % break into sub-matrices for each component
            AA = mat2cell(AA,obj.kConcentration,obj.kConcentration);

            ll = zeros(1,obj.Nc);
            for k=1:obj.Nc
                Lk = obj.Components(k).Concentration.L;
                [~,d] = eig(full(AA{k,k}),full(Lk'*Lk),'qz','vector');
                ll(k) = ng2lambda(d,ng(k));
            end
        end

        function ll = estimateProfileLambda(obj,err,ng)
            
            if nargin < 3 || isempty(ng) % need to estimate ng first
                ng = zeros(1,obj.Nc);
                for j=1:obj.Nc
                    P = obj.Components(j).Profile;
                    switch P.type
                        case 'simple' % no regularization
                            ng(j) = Inf;
                        case 'smooth' % just a little smoothing
                            ng(j) = 0.9*P.maxinfo;
                        case 'realspace' % aggressive smoothing if Ns > 10
                            ng(j) = min([10,0.9*P.maxinfo]);
                        otherwise
                            error('profile type not implemented?');
                    end
                    
                end
                fprintf(1,'estimating profile lambda with ng = %s\n',mat2str(ng));
            end
            
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

            c = obj.concentrations;

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

        function [AA,Ab] = concentrationProblem(obj,I,err)
            % assumes equal errorbars for each q-bin

            w = 1./mean(err,2); % average over q-bins

            nc = numel(obj.Components);
            AA = cell(nc,nc);
            Ab = cell(nc,1);

            Concentrations = [obj.Components.Concentration];
            A = {Concentrations.A};

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
            c = obj.concentrations;
            y = obj.profiles;
            D = I - y(:,notk)*c(:,notk)';
            ck = c(:,k);
            m = ck/(ck'*ck); % coefficients in the sum
            Ik = D*m;
            sigmak = sqrt(err.^2*m.^2);
        end

        function [pk,sigmak] = extractConcentration(obj,I,err,k)
            notk = setdiff(1:obj.Nc,k);
            c = obj.concentrations;
            y = obj.profiles;
            D = I - y(:,notk)*c(:,notk)';
            yk = y(:,k);
            w = 1./mean(err,2); % error weights
            m = (w.^2.*yk)/(yk'*(w.^2.*yk)); % coefficients in the sum
            pk = D'*m;
            sigmak = sqrt(err'.^2*m.^2);
        end

        function H = get.concentrationH(obj)
            B = cell(1,obj.Nc);
            for j=1:obj.Nc
                L = obj.Components(j).Concentration.L;
                B{j} = L*sqrt(obj.concentrationLambda(j));
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

        function p = get.concentrations(obj)
            Concentrations = [obj.Components.Concentration];
            p = zeros(obj.Nx,obj.Nc);
            for j=1:obj.Nc
                p(:,j) = Concentrations(j).A*obj.uConcentration{j};
            end
        end

        function n = get.concentrationNorm(obj)
            n = zeros(1,obj.Nc);
            for j=1:obj.Nc
                n(j) = obj.Components(j).Concentration.norm(obj.uConcentration{j});
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
dd = dd(~isinf(log10(dd))); % fix numberical bug: removes values very close to zero

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
