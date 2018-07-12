# Analysis-tools-for-electrophysiological-somatosensory-cortex-databank

## About
Analysis code from <br />
*da Silva Lantyer, A., Calcini, N., Bijlsma, A., Zeldenrust, F., Scheenen, W. J. J., Celikel, T. (under review) A databank for intracellular electrophysiological mapping of the adult somatosensory cortex. GigaScience* <br />
Please cite this reference when using the frozen noise data.

## Other analysis -- please add description
For instance analysis of VC sawtooth

### About
What does the code do?

### Use
How do I use it?

## Frozen Noise
Using the method as described in: <br />
*Zeldenrust, F., de Knecht, S., Wadman, W. J., Denève, S., Gutkin, B., Knecht, S. De, Denève, S. (2017).  Estimating the Information Extracted by a Single Spiking Neuron from a Continuous Input Time Series.  Frontiers in Computational Neuroscience, 11, 49. doi:10.3389/FNCOM.2017.00049* <br />

Please cite both references when using the frozen noise protocol or data.<br />
See also: https://github.com/fleurzeldenrust/In-vitro-method-for-information-calculation

### About
The analysis code loads a data file from the GigaScience database above, and calculates and plots the mutual information between the hidden state and the spike train. <br />
The generate_input_current code generates a frozen noise input current that can be used in experiments. Note that the hidden_state and input_theory vectors are needed for the information calculation after the experiments!

### Use
see: Frozen Noise / Analysis / Example_analysis_experiments <br /> 
see: Frozen Noise / generate_input_current / Example_generate_input <br /> 
