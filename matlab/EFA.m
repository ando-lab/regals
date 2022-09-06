classdef EFA
    %EFA - Evolving Factor Analysis
    
    properties
        I
        sigma = []
        svd_mode(1,:) char {mustBeMember(svd_mode,{'exact','fast'})} = 'exact'
    end
    properties(Dependent = true)
        Nr
        Nc
    end
    
    methods
        function obj = EFA(I,sigma,varargin)
            % EFA
            
            % assign input arguments
            obj.I = I;
            obj.sigma = sigma;
            
            % assign name, value pairs
            for j=1:2:(numel(varargin)-1)
                obj.(varargin{j}) = varargin{j+1};
            end
        end
        
        function val = get.Nr(obj)
            val = size(obj.I,1);
        end
        
        function val = get.Nc(obj)
            val = size(obj.I,2);
        end
        
        function sv = evolvingFactors(obj,k,direction,skip)
            %evolving_factors
            
            Ncols = obj.Nc;
            
            if nargin < 3 || isempty(direction)
                direction = 'forward';
            end
            
            if nargin < 4 || isempty(skip)
                skip = 1;
            end
            
            sv = NaN*ones(k,Ncols);
            
            for j=1:skip:Ncols
                switch lower(direction)
                    case 'forward'
                        cols = 1:j;
                    case {'reverse'}
                        cols = j:Ncols;
                end
                sv(:,j) = obj.svd(k,cols);
            end
            
        end
        
        function sv = svd(obj,k,cols,normalize,subtract)
            
            if nargin < 2 || isempty(k)
                k = obj.Nc;
            end
            if nargin < 3 || isempty(cols)
                cols = 1:obj.Nc;
            end
            if nargin < 4 || isempty(normalize)
                normalize = true;
            end
            if nargin < 5 || isempty(subtract)
                subtract = true;
            end
            
            sv = NaN*ones(k,1);
            
            % prepare matrix
            A = obj.I(:,cols);
            
            % apply sigma weight
            if ~isempty(obj.sigma) && size(obj.sigma,1) == obj.Nr
                if size(obj.sigma,2) == 1
                    A = A./obj.sigma;
                else
                    A = A./mean(obj.sigma(:,cols),2);
                end
            end
            
            if normalize
                A = A/sqrt(obj.Nr);
            end
            
            % calculate svd
            switch obj.svd_mode
                case 'exact'
                    [~,s,~] = svd(A,'econ');
                    s = diag(s);
                    s = s(1:min(k,numel(s)));
                case 'fast'
                    s = svds(A,k,'largest','Tolerance',1);
            end
            
            % subtract maximum s.v. of Gaussian random matrix
            if subtract && normalize
               s = s - sqrt(numel(cols)/obj.Nr);
            elseif subtract % && ~normalize
               s = s - sqrt(numel(cols));
            end
            
            sv(1:numel(s)) = s;
            
        end
        
        function [y,c,R] = quickRotate(obj,xstart,xend)
            ncomp = numel(xstart);
            
            % do SVD
            w = 1./mean(obj.sigma,2);
            [u,s,v] = svd(obj.I.*w,'econ');
            
            % get the important components
            u = u(:,1:ncomp);
            s = s(1:ncomp,1:ncomp);
            v = v(:,1:ncomp);
            
            % calculate optimal rotation
            R = zeros(ncomp,ncomp);
            
            for n=1:ncomp
                m = false(size(v,1),1);
                m(ceil(xstart(n)):floor(xend(n))) = true;

                vIn = v(m,:); % concentration basis inside peak
                vOut = v(not(m),:); % concentration basis outside peak
                
                % set up least squares problem
                A = vOut;
                B = mean(vIn,1);
                AA = A'*A;
                BB = B'*B;
                
                % large Lagrange multiplier enforces peak normalization
                lambda = 1E6*trace(AA)/trace(BB); 
                
                % solve the normal equations
                R(:,n) = (AA + lambda*(BB))\(lambda*B'*1);
            end

            % calculate I = y*c
            y = (u*s*pinv(R'))./w;
            c = R'*v';
        end
    end
    
    methods(Static = true)
        
        function [xinfl,xc,slope] = fitInflection(s,direction,threshold,window)
            N = size(s,1);
            xinfl = NaN*ones(N,1);
            xc = NaN*ones(N,1);
            slope = NaN*ones(N,1);
            x = 1:size(s,2);
            for j=1:N
                v = s(j,:);
                isIncl = ~isnan(v);
                xj = x(isIncl);
                yj = v(isIncl);
                switch lower(direction)
                    case 'forward'
                        x0 = xj(find(yj>threshold,1,'first'));
                    case 'reverse'
                        x0 = xj(find(yj>threshold,1,'last'));
                    otherwise
                        error('unexpected direction');
                end
                if isempty(x0)
                    continue;
                end
                xc(j) = x0;
                isWindow = (x0 - window/2) <= xj & (x0 + window/2) >= xj;
                npts = nnz(isWindow);
                if npts < 2
                    continue;
                end
                xw = xj(isWindow);
                yw = yj(isWindow);
                A = ones(npts,2);
                A(:,1) = xw(:) - x0;
                c = A\yw(:);
                xinfl(j) = x0 + (1-c(2))/c(1);
                slope(j) = c(1);
            end
            xinfl(xinfl<x(1)) = x(1);
            xinfl(xinfl>x(end)) = x(end);
            
        end
    end
end

