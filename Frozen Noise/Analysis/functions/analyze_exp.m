function Output = analyze_exp(ron, roff, x, input_theory, dt, spiketrain)
% Analyze data from MI experiment
% Zeldenrust, F., de Knecht, S., Wadman, W. J., Denève, S., Gutkin, B., Knecht, S. De, Denève, S. (2017). 
% Estimating the Information Extracted by a Single Spiking Neuron from a Continuous Input Time Series. 
% Frontiers in Computational Neuroscience, 11(June), 49. doi:10.3389/FNCOM.2017.00049
% Please cite this reference when using this method.

% INPUT
% * ron, roff (kHz): switching speed of the hidden state
% * x: array with hidden state values over time
% * input_theory: array with unscaled input current values (output from ANN)
% * dt: binsize recordings (ms)
% * spiketrain: array (same size as x and input_theory) of 0 (no spike) and 1 (spike)

% OUTPUT
% Output-struct with fields:
% * MI_i        : mutual information between hidden state and input current
% * xhat_i      : array with hidden state estimate based on input current 
% * MSE_i       : mean-squared error between hidden state and hidden state estimate based on input current 
% * MI          : mutual information between hidden state and spike train
% * qon, qoff   : spike frequency during on, off state in spike train 
% * xhatspikes  : array with hidden state estimate based on spike train
% * MI          : mean-squared error between hidden state and hidden state estimate based on spike train




Output = struct();

%% Input
[~, ~, Output.MI_i, L_i] = calc_MI_input( ron, roff, input_theory, 0, x, dt );
Output.xhat_i = 1./(1+exp(-L_i));
Output.MSE_i = sum((x - Output.xhat_i).^2);

%% Output
[~, ~, Output.MI, L, Output.qon, Output.qoff] = calc_MI_ideal( ron, roff, spiketrain, x, dt );
Output.xhatspikes = 1./(1+exp(-L));
Output.MSE = sum((x - Output.xhatspikes).^2);


end

%% Helper Functions
    
function [Hxx, Hxy, MI, L, qon, qoff] = calc_MI_ideal( ron, roff, spiketrain, x, dt )
% Calculate the mutual information between hidden state x and output spike
% train (same size vector) assuming a ideal observer that knows ron and roff

% NB Note that if dt in ms, then ron and roff in kHz
% NB Note that information is calculated in bits. For nats use log instead
% of log2

%% Calculate qon and qoff -> w and theta
[spikesup, spikesdown] = reorder_x(x, spiketrain);
spikesup = squeeze(spikesup);
spikesdown = squeeze(spikesdown);

nspikesup = nansum(nansum(spikesup));
qon = nspikesup/(sum(x)*dt);
nspikesdown = nansum(nansum(spikesdown));
if nspikesdown==0
    disp('no down spikes, invent one')
    nspikesdown=1;
end

qoff = nspikesdown/((length(x) - sum(x))*dt);
w = log(qon/qoff);
theta = qon-qoff;

disp(['w = ', num2str(w), '; theta = ', num2str(theta)])


%% Integrate L
input = spiketrain/dt;
L = NaN*ones(size(x));
L(1) = log(ron/roff);
for nn=1:length(x)-1
    L(nn+1) = L(nn)+dLdt_spikes(L(nn), ron, roff, input(nn), w, theta)*dt;
    if abs(L(nn+1))>1000
        disp('L diverges, weigths too large')
        break
    end
end

%% Calculate MI
[Hxx, Hxy, MI] = MIEst(L,x);


end

function [Hxx, Hxy, MI, L] = calc_MI_input( ron, roff, input, theta, x, dt )
% Calculate the mutual information between hidden state x and generated
% input train (same size vector) assuming a ideal observer that knows ron
% and roff and theta.

% NB Note that if dt in ms, then ron and roff in kHz
% NB Note that information is calculated in bits. For nats use log instead
% of log2

%% Integrate L
L = NaN*ones(size(x));
L(1) = log(ron/roff);
for nn=1:length(x)-1
    L(nn+1) = L(nn)+dLdt_input(L(nn), ron, roff, input(nn), theta)*dt;
    if abs(L(nn+1))>1000
        disp('L diverges, weigths too large')
        break
    end
end

%% Calculate MI
[Hxx, Hxy, MI] = MIEst(L,x);

end

function y = dLdt_spikes(L, ron, roff, input, w, theta)
    y = ron*(1+exp(-L))-roff*(1+exp(L)) + w*input - theta;
end

function y = dLdt_input(L, ron, roff, input, theta)
    y = ron*(1+exp(-L))-roff*(1+exp(L)) + input - theta;
end

function [Hxx, Hxy, MI] = MIEst(L, x)
    Hxx = - mean(x)*log2(mean(x)) - (1-mean(x))*log2(1-mean(x));
    Hxy = - mean(x.*log2(sig(L))+(1-x).*log2(1-sig(L)));
    MI =  Hxx - Hxy;
end

function y = sig(x)
    y=1./(1+exp(-x));
end

function [revecsup, revecsdown] = reorder_x(x, ordervecs)
% Reorder the vectors in orderdervec (nvec x length) to x=1 (up) and x=0 (down). 



[nov, Nt] = size(ordervecs);
if nov>Nt
    s = input('Number of vectors larger than number of time steps; transpose? (y/n)','s');
    if strcmp(s, 'y')
        ordervecs = ordervecs';
        [nov, Nt] = size(ordervecs);
    end
end
[novx, ~] = size(x);
if ~(novx == 1)
    x = x';
end

xt1 = [x(1) x];
xt2 = [x x(end)];
xj = xt2-xt1; % this is equal to xj(n) = x(n)-x(n-1)
njumpup = length(find(xj==1));
njumpdown = length(find(xj==-1));

%% reorder
if njumpup>0 && njumpdown>0
    firstjump = find(abs(xj)==1);
    firstjump = firstjump(1);
    revecsup = NaN*ones(nov,njumpup+1 ,round(10*Nt/njumpup));
    revecsdown = NaN*ones(nov,njumpdown+1 ,round(10*Nt/njumpdown));
    [~, ~, size3] = size(revecsdown);

    if x(firstjump) == 1
        up = 1;
        down = 0;
        revecsup(:,1,1) = ordervecs(:,firstjump);
    elseif x(firstjump)==0
        up = 0;
        down = 1;
        revecsdown(:,1,1) = ordervecs(:,firstjump);
    else
        error('first jump not properly defined')
    end
    tt = 1;
    tmaxup = 1;
    tmaxdown = 1;

    for nn=firstjump+1:Nt
        try
            jump = x(nn)-x(nn-1);
        catch err
            disp('size ordervecs not the same as size x!')
            keyboard        
        end
        if jump == 0
            tt = tt+1;
            if x(nn)==1
                % up state
                if tt>tmaxup
                    tmaxup = tt;
                end
                revecsup(:,up,tt) = ordervecs(:,nn);
            elseif x(nn)==0
                % down state
                if tt>tmaxdown
                    tmaxdown = tt;
                end
                revecsdown(:,down,tt) = ordervecs(:,nn);
            else 
                error('something went wrong: x not 0 or 1')
            end
        elseif jump == 1
            % jump up
            tt = 1;
            up = up+1;
            if x(nn)==1
                revecsup(:,up,tt) = ordervecs(:,nn);
            else
                error('something went wrong: jump up but x is not 1')
            end
        elseif jump == -1
            % jump down
            tt = 1;
            down = down+1;
            if x(nn)==0
                revecsdown(:,down,tt) = ordervecs(:,nn);
            else
                error('something went wrong: jump down but x is not 0')
            end
        else
            error('something went wrong: no jump up or down')
        end
        if tt>size3-1
            % tt will run out of matrix size
            error('choose larger starting matrix')
        end
    end


    revecsup = revecsup(:,1:up, 1:tmaxup);
    revecsdown = revecsdown(:,1:down,1:tmaxdown);
else
    if njumpup<1
        disp('no jumps up; reordering not possible')
        revecsup = NaN;
        revecsdown = NaN;
    end
    if njumpdown<1
        disp('no jumps down; reordering not possible')
        revecsup = NaN;
        revecsdown = NaN;
    end
end

end
