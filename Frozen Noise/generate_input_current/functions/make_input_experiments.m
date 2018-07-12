function [input_current, input_theory, hidden_state] = make_input_experiments(baseline, amplitude_scaling, tau, factor_ron_roff, mean_firing_rate, sampling_rate, duration, seed)
% Make input current based on a artificial network responding to a hidden state
% The method is described in the following paper:
% Zeldenrust, F., de Knecht, S., Wadman, W. J., Denève, S., Gutkin, B., Knecht, S. De, Denève, S. (2017). 
% Estimating the Information Extracted by a Single Spiking Neuron from a Continuous Input Time Series. 
% Frontiers in Computational Neuroscience, 11(June), 49. doi:10.3389/FNCOM.2017.00049
% Please cite this reference when using this method.


% INPUT
% * baseline, amplitude_scaling: scaling of the stimulus (for instance the current needed to keep the
% membrane potential at a fixed value); note that this scales the 'arbitrary values' of the
% Bayesian network to for instance pA, so it is a bit arbitrary
% * tau (ms): switching speed of the hidden state
% * mean_firing_rate(kHz): mean firing rate of the artificial neurons
% * sampling rate (kHz): of the experimental setup (injected current) 
% * duration (ms): length of the experiment
% * seed (optional) for seeding the random number generator
% NB Make sure that you save the hidden state with the experiments, it is
% essential for the information calculation!

% OUTPUT: 
% * input current: (1xNtime) array with current values
% * input_theory: unscaled, theoretical input
% * hidden state: (1xNtime) array with hidden state values

% Example use:
% [input_current, hidden_state] = make_input_experiments(0, 1000, 20/3000,40/3000, (0.5)/1000, 1, 20000);
% This creates an input current with baseline 0 pA, amplitude 1000 pA, tau=50 ms, the mean firing rate of neurons in the artificial network is 0.5 Hz, sampling rate of 1 kHz, 20000 ms long (you will need at least about 20 s for a good estimate. 

addpath('Classes')

%% Parameters
if nargin == 7
    rng('shuffle');
    seed = rand;
end

dt = 1./sampling_rate;

%% Fixed Parameters
N = 1000;                   % # neurons in artificial network (interchangable with meanq)
tau_exponential_kernel = 5; % ms
alpha = sqrt(1./8.);        % relation mean and std of the firing rates of the artificial neurons (firing rates should be positive)


stdq = alpha*mean_firing_rate;

%% Create input from artifical network
ron = 1./(tau*(1+factor_ron_roff));
roff = factor_ron_roff*ron;
input_bayes             = Input;
input_bayes.fHandle     = @Input.markov;
input_bayes.dt          = dt;
input_bayes.T           = duration;
input_bayes.kernel      = 'exponential';    % EPSC shape
input_bayes.kerneltau   = tau_exponential_kernel;    % in ms
input_bayes.ron         = ron;       % in kHz,  on-switch rate
input_bayes.roff        = roff;
input_bayes.seed        = seed;
input_bayes.xseed       = seed;
[input_bayes.qon, input_bayes.qoff] = input_bayes.create_qonqoff_balanced(input_bayes, N, mean_firing_rate, stdq, seed);
input_bayes             = input_bayes.generate;

%% Scale for experiments
input_theory = input_bayes.input;
input_current = amplitude_scaling*input_bayes.input+baseline;
hidden_state = input_bayes.x;
