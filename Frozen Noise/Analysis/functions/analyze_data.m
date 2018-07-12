function Data = analyze_data(Data, Params)

try
    windowtype = Params.MIanalysis.windowtype;
    windowsize = Params.MIanalysis.windowsize;
catch
    if ~(strcmp(windowtype, 'all') || strcmp(windowtype, 'dependsontau'))
        error('Please give a windowsize')
    end
end
settings = Data.settings;


dt = 1/settings.sampling_rate;
input_current = Data.input_current;
hidden_state = Data.hidden_state;
membrane_potential = Data.membrane_potential;
    
if strcmp(windowtype, 'all')
    indexwindow = length(hidden_state);
elseif strcmp(windowtype, 'timewindow')
    indexwindow = round(Params.MIanalysis.windowsize/dt);
elseif strcmp(windowtype, 'indexwindow')
    indexwindow = Params.MIanalysis.windowsize;
elseif strcmp(windowtype, 'dependsontau')
    if ~isfield(Params.MIanalysis, 'windowsize')
        windowsize50 = 20000; % 20 s window for analysis at tau = 50ms
    else
        windowsize50 = Params.MIanalysis.windowsize;
    end
    if Data.settings.tau <= 50
        indexwindow = round(windowsize50/dt);
    else
        wfac = Data.settings.tau/50;
        indexwindow = round(wfac*round(windowsize50/dt));
    end
else
    error('Please give an appropriate window type')
end
    
%% Fixed Parameters

ron = 1./(settings.tau*(1+Params.MIanalysis.factor_ron_roff));
roff = Params.MIanalysis.factor_ron_roff*ron;
input_theory = (input_current - (settings.baseline))./ (settings.amplitude_scaling);

Ntime = length(Data.membrane_potential);
Nwindow = floor(Ntime/indexwindow);
    
%% Find spikes
searchthreshold = 0; %mV
[~, Data.spikeindices] = findspikes(dt, membrane_potential, searchthreshold, 0);
spiketrain = zeros(size(hidden_state));
spiketrain(Data.spikeindices) = 1;
    
    
                
for nw = 1:Nwindow
    window_temp = (nw-1)*indexwindow+1:nw*indexwindow;

    hidden_state_temp       = hidden_state(window_temp);
    input_theory_temp       = input_theory(window_temp);
    spiketrain_temp         = spiketrain(window_temp);

    Data.firing_rate{1, nw} = sum(spiketrain_temp)/(indexwindow*dt/1000);

    Data.Analysis{1, nw} = analyze_exp(ron, roff, hidden_state_temp', input_theory_temp, dt, spiketrain_temp');
    Data.Analysis{1, nw}.FI = Data.Analysis{1,nw}.MI/Data.Analysis{1,nw}.MI_i;
    if abs(Data.Analysis{1, nw}.FI)>2
        disp('FI>2, should not be possible')
    end
    if Data.Analysis{1, nw}.MI_i<0
        disp('MI in input negative')    
    end
end


            
end
