# Warm - UP
Little set of code lines (warm-up) that helps you passing the 1st technical interview : ) Each topic is only a few lines of code
- default and optional arguments
- class methods
- super
- lambda
- map
- recursion
- multiprocessing
- decorators



# Simple GUI 
How do build a simple GUI (beginner level). PDF from ppt attached. GUI "My_GUI.py" calls the function "my_data_analysis.py" which uses "poissfit.py" for analyzing peptide count data (very basic!). Data set peptides.xlsx was provided with permisssion from Dr Muriel Teeseling. The functions "my_data_analysis.py" and "poissfit.py" are called from "MyModule.py"



# Bayesian Signal Detection 
Code for detecting a signal from count data (photon time of arrival, ToA) as an alternative to FFT (which only works well for actual light curves). See Gregory & Loredo 1992. This algorithm has been used for X-Ray pulsar timing, see the papers of Hambaryan, Hohle et al., 201X <br/>
run eg:<br/>

from BayesianSignalDetect import * <br/>
T = CreateSignal(N = 5000, w = 0.1, noiseratio = 0.5)                     #creates ToA (T) and creates plot of binned light curve <br/>
FFT(T)                                                                    #FFT barely, if at all, finds w. <br/> 
S = SignalDetect(T)<br/> 
[Omega, P] = S.FindFrequency()                                            #should give you a nice periodogram after 40s <br/>											
See also "Results.pdf".

Note 1: code uses multiprocessing, i. e. performance depends on n_cpu <br/>
Note 2: photon count is always positive (!), i. e. a sin wave of w = 1 will be detected as w = 2 with another peak at w = 1 <br/>