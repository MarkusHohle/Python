# Bayesian Signal Detection 
Code for detecting a signal from count data (photon time of arrival, ToA) as an alternative to FFT (which only works well for actual light curves). See Gregory & Loredo 1992. This algorithm has been used for X-Ray pulsar timing, see the papers of Hambaryan, Hohle et al., 201X <br/>
run eg:<br/>

from BayesianSignalDetect import * <br/>
T = CreateSignal(N = 5000, w = 0.44, noiseratio = 0.1)                    #creates ToA (T) and creates plot of binned light curve <br/>
FFT(T)                                                                    #FFT barely, if at all, finds w. <br/> 
S = SignalDetect(T)<br/> 
[Omega, P] = S.FindFrequency()                                            #should give you a nice periodogram after 40s <br/>											
See also "Results.pdf".

Note 1: code uses multiprocessing, i. e. performance depends on n_cpu <br/>
Note 2: photon count is always positive (!), i. e. a sin wave of w = 1 will be detected as w = 2 with another peak at w = 1 <br/>

Default and keyword arguments: <br/>

S = SignalDetect(T, Range_phi = [0, 2*np.pi], dphi = 0.01, MaxM = 20, **w_start = 0.1, **w_end = 0.5, **dw = 0.001, **dt = 1) <br/>

T	   :   Signal <br/>
Range_phi  :   phase <br/>
dphi       :   increment for phase<br/>
MaxM       :   max number of bins for light curve<br/>
#Opts:<br/>
w_start    :   start value for frequency in grid search<br/>
w_end      :   stop value for frequency in grid search<br/>
dw         :   increment for (angular) frequency<br/>
dt         :   time resolution of detector






