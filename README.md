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

`
MVals - array containing measurement values
`

`
SVals - array containing sparsity values
`

The above scripts will generate the plots as shown in the paper with each block showing the normalized error between the true and recovered signal. If you find this work useful, feel free to use the following citation for the same -

` ` `
@INPROCEEDINGS{9413470,  

author={Vora, Jian and Rajwade, Ajit},  
booktitle={ICASSP 2021 - 2021 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},   
title={Compressive Signal Recovery Under Sensing Matrix Errors Combined With Unknown Measurement Gains},   
year={2021}, 
pages={5105-5109},  
doi={10.1109/ICASSP39728.2021.9413470}}
` ` `
