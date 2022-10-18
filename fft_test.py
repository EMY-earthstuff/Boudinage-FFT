#code snippet to generate different FFT plots of
#mapped contacts extracted for analysis. Note: drops
#synthetic and test contacts from input files.


import numpy as np
from scipy import interpolate
from scipy.fft import fft, ifft, fftshift, fftfreq
from scipy import signal
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy.signal import find_peaks

#change id to load correct file
Run_id = 'normal'


#taper function
def taper(signal):
    xs = np.linspace(0, len(signal), len(signal)+1)
    st = signal[0]
    nd = signal[-1]
    di = nd - st
    m = (signal[-1] - signal[-2])
    b_pos =  signal[-1] - m*xs[-1]
    b_neg = signal[-1] + m*xs[-1]
   
    X = xs[-1]
    Y = signal[-1]
    out = np.zeros(signal.shape)
    out = out + signal
   
    if nd < st:
        while Y < st:
            X+=1
            if m > 0:
                Y = m*X + b_pos
            else:
                Y = -m*X + b_neg
            out = np.append(out, Y)
    else:
        while Y > st:
            X+=1
            if m < 0:
                Y = m*X + b_pos
            else:
                Y = -m*X + b_neg
            out = np.append(out, Y)
           
    return out
   
#scaling function    
def scaling_A_2_lambda(a, f):
    peaks = find_peaks(np.abs(a))
    highA = np.abs(a)[peaks[0][0]]
    highf = 1/f[peaks[0][0]]
   
    m = highA/highf
   
    return m
   
###plotting loops    
   
#bunch of empty lists to dump values in once modeified by functions
sums = []
sums_r = []
sums_a = []
maxs_total = []
maxs_reduced = []

#gen from text the files
boudins = np.genfromtxt('boudin_' + Run_id + '.txt', skip_header=1, delimiter=',')

#some funky code to plot only refs as thick visible lines
normal_refs = [0,1,3]
ptach_refs = [0,1,2]


normal_cs = ['darkred', 'indianred', 'orangered', 'darkorange', 'orange', 'goldenrod', 'gold', 'magenta']
ptach_cs = ['SlateBlue2', 'RoyalBlue1', 'blue2', 'SteelBlue1', 'DeepSkyBlue2', 'SkyBlue2', 'SlateGray3', 'SlateGray4']

#set conditions for plotting right refe contacts
if Run_id == 'ptach':
    cs = ptach_cs
    rangeval = 6
    refs = ptach_refs
else:
    cs = normal_cs
    rangeval = 8
    refs = normal_refs

plt.figure()

for i in range(0, rangeval):
   

    #get index locations of each test boudin and ids
    ptach = np.where(boudins[:,4] == i)
    #get length
    length = boudins[ptach][0,3]
    #get regular intervals and interpolation to get regular points along contact
    samples = 500 #sample number

    ptx = np.linspace(boudins[ptach][:,5].min(), boudins[ptach][:,5].max(), samples)
    
    ptws = ptx[-1] - ptx[0] 

    pt = interpolate.interp1d(boudins[ptach][:,5], boudins[ptach][:,6])

    #detrend & taper & and fetch signal spacing
    ptdt = signal.detrend(pt(ptx))
    pttp = taper(ptdt)



    #fft, adapted from https://pythontic.com/visualization/signals/fouriertransform_fft
    #normalize amplitude
    ptfft = np.fft.fft(pttp)/len(pttp)
    #exclude sampling frequency    
    ptfft = ptfft[range(int(len(pttp)/2))]
    #get frequencies
    tpCount     = len(pttp)
    values      = np.arange(int(tpCount/2))
    timePeriod  = tpCount/samples
    frequencies = values/timePeriod
    
    
    ##plots
    
    if  refs.count(i) == 1:
        apha = 1.
        liw = 3.
    else:
        apha = 0.25
        liw = 1.
   
    for j in range(0, len(ptfft[0:len(frequencies)])):
        if 1/frequencies[j] > 2*ptws:
            frequencies[j] = np.nan
            
    #frequency vs amplitudes
    plt.plot(1/frequencies, np.abs(ptfft), '-', alpha = apha, lw = liw)

    plt.title(str(i))
    plt.xscale('log')

    plt.legend()
    plt.grid('True')    
    
    