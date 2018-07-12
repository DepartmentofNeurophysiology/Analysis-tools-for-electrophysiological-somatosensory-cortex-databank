addpath('functions')


%% Parameters
baseline = 0;           % pA
amplitude = 700;        % pA
tau = 50;               % ms
factor_ron_roff = 2;    % how much more often the hidden state is in the off-state than in the on-state
mean_firing_rate = (0.5)/1000; % firing rate of neurons in the artificial neural net
sampling_rate = 5;      % kHz
dt = 1/sampling_rate; 
duration = 20000;       % ms

%% Make
[input_current, input_theory, hidden_state] = make_input_experiments(baseline, amplitude, tau, factor_ron_roff, mean_firing_rate, sampling_rate, duration);
