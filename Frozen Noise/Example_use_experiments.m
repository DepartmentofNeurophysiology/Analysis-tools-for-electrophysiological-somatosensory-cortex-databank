filename = 'example_data.mat';  % put the full path to a frozen noise file here.

Settings.MIanalysis.windowtype = 'timewindow';
Settings.MIanalysis.windowsize = 20000;          % ms, time window for analysis; use 20000 for tau = 50 ms and 100000 for tau = 250 ms;   
Settings.MIanalysis.factor_ron_roff = 2;         % how much more often in off state than on state

Data = analyze_plot(filename, Settings);
