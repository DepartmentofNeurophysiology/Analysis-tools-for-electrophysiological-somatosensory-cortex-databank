# Analysis-tools-for-electrophysiological-somatosensory-cortex-databank

## About
Analysis code for the frozen noise data from <br />
*da Silva Lantyer, A., Calcini, N., Bijlsma, A., Zeldenrust, F., Scheenen, W. J. J., Celikel, T. (2018) A databank for intracellular electrophysiological mapping of the adult somatosensory cortex. GigaScience* <br />

## Frozen Noise
Using the method as described in: <br />
*Zeldenrust, F., de Knecht, S., Wadman, W. J., Denève, S., Gutkin, B., Knecht, S. De, Denève, S. (2017).  Estimating the Information Extracted by a Single Spiking Neuron from a Continuous Input Time Series.  Frontiers in Computational Neuroscience, 11, 49. doi:10.3389/FNCOM.2017.00049* <br />

Please cite both references when using the frozen noise protocol or data.<br />
See also: https://github.com/fleurzeldenrust/In-vitro-method-for-information-calculation

### About
In the Analysis folder, the function 'analyze_plot' loads a data file from the GigaScience database above, and calculates and plots the mutual information between the hidden state and the spike train. An example script is given: 'Example_analysis_experiments'.<br />
In the generate_input_current folder, it is shown how to generate a frozen noise input current that can be used in experiments. Note that the hidden_state and input_theory vectors are needed for the information calculation after the experiments! In an example script 'Example_generate_input' it is shown how this code is used.

### Use
see: Frozen Noise / Analysis / Example_analysis_experiments <br /> 
see: Frozen Noise / generate_input_current / Example_generate_input <br /> 
