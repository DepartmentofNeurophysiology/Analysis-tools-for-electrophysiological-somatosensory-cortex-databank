function  Data = analyze_plot(filename, Settings)


%% Initialize 
f = filesep;

if nargin < 2
    % Analysis
    Settings.MIanalysis.windowtype = 'dependsontau';
    Settings.MIanalysis.windowsize = 20000;             % ms, time window for analysis of tau=50 ms files, other tau scaled accordingly.    
    Settings.MIanalysis.factor_ron_roff = 2;            % how much more often in off state than on state    
end

addpath(genpath(['.' f 'functions']));

%% Get data
disp(['Loading ' filename])
Data = load(filename);
Data.membrane_potential = Data.cell_response*1000;
Data.input_current = Data.input_current*1e12;
Data = rmfield(Data,'cell_response');


%% Analyze
disp('Analysing')
Data = analyze_data(Data, Settings);

%% Plot
disp('Plotting')
dt = 1/Data.settings.sampling_rate;
time = dt*(1:length(Data.hidden_state));
timewindow = [3000, 5000];
nn = round(timewindow(1)/dt):round(timewindow(2)/dt);


figure
legendtext = {};
        
color = [];
symbol = [];
try
    if (strcmp(Data.settings.condition, 'aCSF') || strcmp(Data.settings.condition,'aCSf'))
        color = 'k';
        symbol = '+';
        legendtext(end+1) = {Data.settings.condition};
    elseif strcmp(Data.settings.condition, 'Dop')
        color = 'r';
        symbol = 'd';
        legendtext(end+1) = {'Dopamine'};
    elseif strcmp(Data.settings.condition, 'D1') 
        color = 'b';
        symbol = '*';
        legendtext(end+1) = {'D1 antagonist'};
    elseif strcmp(Data.settings.condition, 'D2') 
        color = 'g';
        symbol = 'o';
        legendtext(end+1) = {'D2 antagonist'};
    elseif strcmp(Data.settings.condition, 'D1ago')
        color = 'c';
        symbol = 'x';
        legendtext(end+1) = {'D1 agonist'};
    elseif strcmp(Data.settings.condition, 'D2ago')
        color = 'm';
        symbol = 's';
        legendtext(end+1) = {'D2 agonist'};
    elseif strcmp(Data.settings.condition, 'BIC')
        color = 'y';
        symbol = '^';
        legendtext(end+1) = {'Bicuculline'};
    else
        color = 'k';
        symbol = '.';
        legendtext(end+1) = {'Undefined'};
    end
catch
    disp('Conditions not properly defined')
    color = 'k';
    symbol = '.';
    legendtext(end+1) = {'Undefined'};
end


subplot(3,1,1)
hold all
plot((time(nn)-time(nn(1))), Data.input_current(nn), 'LineWidth',2)
yl = get(gca, 'YLim');
fcvec = [.8 .8 .8];

jumps = [ Data.hidden_state; nan]-[nan; Data.hidden_state];
jumpup = find(jumps == 1);
jumpdown = find(jumps == -1);
jumpup = jumpup(jumpup>nn(1));
jumpup = jumpup(jumpup<nn(end));
jumpdown = jumpdown(jumpdown>nn(1));
jumpdown = jumpdown(jumpdown<nn(end));

if jumpup(1)>jumpdown(1)
    % rectangle at the edge beginning
    x1 = 0;
    x2 = dt*jumpdown(1)-time(nn(1));
    p=patch([x1 x2 x2 x1],[yl(1) yl(1) yl(2) yl(2)],fcvec);
    set(p,'FaceAlpha',0.5, 'EdgeColor','none');
    jumpdown(1) = [];
end
if jumpup(end)>jumpdown(end)
    % recrangle at the edge end
    x1 = dt*jumpup(end)-time(nn(1));
    x2 = time(nn(end))-time(nn(1));
    p=patch([x1 x2 x2 x1],[yl(1) yl(1) yl(2) yl(2)],fcvec);
    set(p,'FaceAlpha',0.5, 'EdgeColor','none');     
    jumpup(end) = [];
end
for jnr = 1:length(jumpup)
    x1 = dt*jumpup(jnr)-time(nn(1));
    x2 = dt*jumpdown(jnr)-time(nn(1));
    p=patch([x1 x2 x2 x1],[yl(1) yl(1) yl(2) yl(2)],fcvec);
    set(p,'FaceAlpha',0.5, 'EdgeColor','none');
end
ylim(yl)
ylabel('current (pA)')
title('Input current')
xlabel('time (ms)')
set(gca, 'FontSize',16)
grid on
box on


subplot(3,1,2)
hold all
plot(time(nn)-time(nn(1)), Data.membrane_potential(nn), color, 'LineWidth',2)
if isfield(Data, 'spikeindices')
    plot(dt*(Data.spikeindices-nn(1)), Data.membrane_potential(Data.spikeindices), ['o' color])
end
xlim([timewindow(1)-time(nn(1)), timewindow(2)-time(nn(1))])
xlabel('time (ms)')
ylabel('V_m (mV)')
title('Membrane potential')
set(gca, 'FontSize',16)
grid on
box on

%% Information plots
if isfield(Data, 'firing_rate')
    fr = cell2mat(Data.firing_rate);
    Nwindow = length(fr);
    FI = nan(1,Nwindow);
    MI_i = nan(1,Nwindow);
    MI_o = nan(1,Nwindow);
    if isfield(Data, 'Analysis')
        for nw = 1:Nwindow
            FI(nw) = Data.Analysis{nw}.FI;
            MI_i(nw) = Data.Analysis{nw}.MI_i;
            MI_o(nw) = Data.Analysis{nw}.MI;
        end
    end

    subplot(3,3,7)
    hold all
    plot(MI_i, MI_o, [color symbol]);
    xlabel('MI_{input} (bits)')
    ylabel('MI_{output} (bits)')
    ylim([-.5 1])
    xlim([-.5 1])
    title(['Input versus output information'])
    set(gca, 'FontSize',16)
    grid on
    box on
    
    subplot(3,3,8)
    hold all
    plot(fr, FI, [color symbol]);
    xlabel('firing rate (Hz)')
    ylabel('FI (MI_o / MI_i, unitless)')
    ylim([-.5 1])
    title(['Fraction of transferred information'])
    set(gca, 'FontSize',16)
    grid on
    box on
    
    subplot(3,3,9)
    hold all
    plot(fr*Data.settings.tau/1000, FI, [color symbol]);
    xlabel('normalized firing rate (r*\tau, unitless)')
    ylabel('FI (MI_o / MI_i, unitless)')
    ylim([-.5 1])
    title(['Fraction of transferred information'])
    set(gca, 'FontSize',16)
    grid on
    box on


end




    

