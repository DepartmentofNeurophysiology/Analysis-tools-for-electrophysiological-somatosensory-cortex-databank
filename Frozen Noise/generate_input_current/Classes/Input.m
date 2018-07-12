classdef Input % value class
    % NB time (dt,T) in ms, freq in Hz, maar qon en qoff in MHz
    % always first 50 ms silent, then 50 ms noise, then rest
    properties
        % for all
        dt          % in ms
        T           % in ms
        fHandle     % type input
        seed        % for seeding the random generator
        input       % final vector
        
        % for Markov
        ron         % in MHz
        roff        % in MHz
        qon         % in MHz
        qoff        % in MHz
        kernel
        kerneltau   % in ms
        xseed
        x
        xfix        % generate input for given x
        
    end
    
    properties (Dependent, Hidden)
        tvec
        length
        tau
        P0
        theta
        w
    end
    
    methods
        function b = Input()
        end
        
        function tv=get.tvec(b)
            tv=b.dt:b.dt:b.T;
        end
        
        function le=get.length(b)
            le=length(b.tvec);
        end
        
        function b = generate(b)
            [b.input, b.x]=b.fHandle(b);
        end
        
        function p = get.tau(b)
            if ~isempty(b.ron)
                if ~isempty(b.roff)
                    p = 1./(b.ron+b.roff);
                else
                    disp('tau not defined for this type of input')
                    p = NaN;
                end
            else
                disp('tau not defined for this type of input')
                p = NaN;
            end
        end
        
        function p = get.P0(b)
            if ~isempty(b.ron)
                if ~isempty(b.roff)
                    p = b.ron./(b.ron+b.roff);
                else
                    disp('P0 not defined for this type of input')
                    p = NaN;
                end
            else
                disp('P0 not defined for this type of input')
                p = NaN;
            end
        end
        
        function p = get.theta(b)
            if ~isempty(b.qon)
                if ~isempty(b.qoff)
                    p = sum(b.qon-b.qoff);
                else
                    disp('theta not defined for this type of input')
                    p = NaN;
                end
            else
                disp('theta not defined for this type of input')
                p = NaN;
            end
        end
        
        function p = get.w(b)
            if ~isempty(b.qon)
                if ~isempty(b.qoff)
                    p = log(b.qon./b.qoff);
                else
                    disp('w not defined for this type of input')
                    p = NaN;
                end
            else
                disp('w not defined for this type of input')
                p = NaN;
            end
        end
        
    end
    
    methods (Static)
        
        function [qon, qoff] = create_qonqoff(~, mutheta, N, alphan, regime, qseed)

            if nargin==5
                rng('shuffle')
            elseif nargin == 6
                rng(qseed)
            else
                disp('not correct number of input arguments')
                qon = nan;
                qoff = nan; 
                return
            end

            qoff = randn(1,N);
            if N>1
                qoff = qoff/std(qoff);
            end
            qoff = qoff - mean(qoff);
            qon = randn(1,N);
            if N>1
                qon = qon/std(qon);
            end
            qon = qon - mean(qon);
            if regime == 1
                % Coincidence regime
                qoff = (alphan*qoff+1)*mutheta/N;
                qon = (alphan*qon+2)*mutheta/N;
            else
                % Push-pull regime
                qoff = (alphan*qoff+1)*mutheta/sqrt(N);
                qon = (alphan*qon+1+1/sqrt(N))*mutheta/sqrt(N);
            end
            qoff(qoff<0) = abs(qoff(qoff<0));
            qon(qon<0) = abs(qon(qon<0));
        end
        
        function [qon, qoff] = create_qonqoff_balanced(~,N,  meanq, stdq, qseed)


            if nargin==4
                rng('shuffle')
            elseif nargin == 5
                rng(qseed)
            else
                disp('not correct number of input arguments')
                qon = nan;
                qoff = nan; 
                return
            end

            qoff = randn(1,N);
            if N>1
                qoff = qoff/std(qoff);
            end
            qoff = qoff - mean(qoff);
            qon = randn(1,N);
            if N>1
                qon = qon/std(qon);
            end
            qon = qon - mean(qon);
            
            qoff = stdq*qoff+meanq;
            qon = stdq*qon+meanq;
            
            qoff(qoff<0) = abs(qoff(qoff<0));
            qon(qon<0) = abs(qon(qon<0));
        end
        
        function [qon, qoff] = create_qonqoff_balanced_uniform(~,N,  minq, maxq, qseed)


            if nargin==4
                rng('shuffle')
            elseif nargin == 5
                rng(qseed)
            else
                disp('not correct number of input arguments')
                qon = nan;
                qoff = nan; 
                return
            end

            qoff = minq + (maxq-minq).*rand(1,N);
            qon = minq + (maxq-minq).*rand(1,N);
        end
        
        function [ip, xs]=markov(b)
            % takes vector with qon and qoff values + values for ron and roff,
            % generates input and x if xfix empty, generates input with xfix otherwise;
            
            if isempty(b.xseed)
                rng('shuffle')
            else
                rng(b.xseed)
            end
            
            ni=length(b.qon);
            nt=length(b.tvec);
                        
            w = log(b.qon./b.qoff);
            
            % generate x
            if isempty(b.xfix)
                P0=b.ron/(b.ron+b.roff);
                xs=zeros(size(b.tvec));
                % initial value
                i=rand;
                if i<P0
                    xs(1)=1;
                else
                    xs(1)=0;
                end

                % make x
                for n=2:length(b.tvec)
                    i=rand;
                    if xs(n-1)==1
                        if i<b.roff*b.dt
                            xs(n)=0;
                        else
                            xs(n)=1;
                        end
                    else
                        if i<b.ron*b.dt
                            xs(n)=1;
                        else
                            xs(n)=0;
                        end
                    end

                end 
            else
                xs = b.xfix;
            end

            % make spike trains (implicit) 
            stsum = zeros(1,nt);
            if ~isempty(b.kernel)
                if strcmp(b.kernel, 'exponential')
                    tfilt=0:b.dt:5*b.kerneltau;
                    kernelf=exp(-tfilt/b.kerneltau);
                    kernelf=kernelf/(b.dt*sum(kernelf));
                elseif strcmp(b.kernel, 'delta')
                    kernelf = 1./b.dt;                    
                end
            end
            xon = find(xs==1);
            xoff = find(xs==0);
            if isempty(b.seed)
                rng('shuffle')
            else
                rng(b.seed)
            end
            for k=1:ni
                
                randon = rand(size(xon));
                randoff = rand(size(xoff));
                sttemp = false(1,nt);
                sttempon = false(size(xon));
                sttempoff = false(size(xoff));
                
                sttempon (randon <b.qon(k) *b.dt) = 1;
                sttempoff(randoff<b.qoff(k)*b.dt) = 1;
                
                sttemp(xon)=sttempon;
                sttemp(xoff)=sttempoff;
                
                stsum = stsum+w(k)*sttemp;
            end   
            if ~isempty(b.kernel)
                stsum=conv(stsum, kernelf, 'full'); 
            end
            stsum = stsum(1:nt);
            ip = stsum;
        end 
        
    end
end



