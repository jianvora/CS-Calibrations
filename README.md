# Compressive Signal Recovery Under Sensing Matrix Errors Combined With Unknown Measurement Gain

Code for our ICASSP 2021 paper which can be found [here](https://ieeexplore.ieee.org/abstract/document/9413470)

The above is a module to blindly calibrate errors (both in the sensing matrix along with multiplicative gains) during compressive signal acquisition. SensorGain.m is the script for canonical sparsity while ObjectPose.m when the signal is sparse in the Haar Wavelet basis.

Some parameters to better understand the above scripts:

`
N - length of the signal to be acquired
`
`
r - maximum frequency perturbation value in MRI acquisition
`
`
noisefrac - fraction of noise in the measurements (for simulation purposes)
`
`
r_gain - maximum gain perturbation
`
`
numdeltas - number of unique frequency perturbation parameters
`
`
numgains - numer of unique gain perturbation parameters
`
