# acquire spectrum from Ocean Optics USB200
# and process data
# 22-May-2024 J.Beale

import numpy as np
from matplotlib import pyplot as plt 
from seabreeze.spectrometers import Spectrometer
import time  # for time-date and real-time delay
import bottleneck as bn  # boxcar filter
import os  # combine dir + filename into path
import pandas as pd  # write .csv output file
from scipy.signal import find_peaks

outDir = "/home/john/Documents/Spectrometer"  # directory for output data files

integrationTime = 50000 # spectrometer integration time in microseconds
averages = 50          # how many spectra to average together
boxcar = 1             # boxcar-filter this many adjacent pixels (always odd)
# xRange = (360, 1000)    # display this range in nm

#yRange = (0, 4200)      # this range in intensity units (12-bit ADC)
yRange = (0, 70)      # this range in intensity units (12-bit ADC)

#  peak-finding parameters
wSize= 11  # rolling average window size
pkSearchWin = 41 # array element range, must be odd
pkStart = 500  # minimum wavelength peak of interest
pkStop = 800   # maximum wavelength peak of interest
xRange = (pkStart, pkStop)    # display this range in nm


# ==================================================================

# replace all NaN values with 0
def rmN(A):
    A[np.isnan(A)] = 0
    return(A)

# boxcar-filter array A with window size w
# w should be an odd number, otherwise result has offset
def filt(A, w):
    b = bn.move_mean(A, window=w) # should left-shift result by (n-1)/2
    b = rmN(b)
    res = np.roll(b, -int((w-1)/2) )
    return(res)

# get a spectrum which is the average of N acquisitions, boxcar Avg pixels together
def getSpec(N, Avg):    
    A = spec.intensities()
    for i in range(N-1):
        A += spec.intensities()
    A = A / N  # convert the sum to an average
    if (Avg > 1): # average together neighboring pixels?
        A = filt(A, Avg)
    return(A)

# display peak location
def printPeaks(pkIdx, A, pos, start, stop):
    dCount = 0
    for i in pkIdx:  # display position of peaks in nm
        pk = pos[i]
        amp = A[i]
        if (pk >= start) and (pk <= stop):
            print("pk = %5.3f  size: %5.3f" % (pk, amp))
            dCount += 1
    print("Peaks: %d" % dCount)
    

# ==================================================================
# main code starts here

spec = Spectrometer.from_first_available()
spec.integration_time_micros(integrationTime)

nm = spec.wavelengths()

# plt.axis([365,1000,0,4096])
plt.ion()
plt.show()

baseline = getSpec(averages, boxcar)  # get baseline reference

while True:


    A = getSpec(averages, boxcar) - baseline
    As = np.sqrt(np.clip(A,0,4095))
    plt.clf()  # clear current plot figure
    plt.plot(nm, As)
    plt.xlim(xRange)
    plt.ylim(yRange)
    plt.title("Spectrum")
    plt.xlabel('wavelength, nm')
    plt.ylabel('intensity')
    plt.grid(axis='both')

    # find peaks ============================
    Afilter = A
    specFF = filt(Afilter, pkSearchWin)
    for i in range(2):
        specFF = filt(specFF, pkSearchWin)

    specHP = Afilter - specFF
    specHPF = filt(specHP, pkSearchWin*2)
    specHPF = filt(specHPF, pkSearchWin*2)
    specHPP = specHP - specHPF  # 2nd stage highpass filter
    specHPP = np.clip(specHPP,0.002,1000)  # clip all negative values to zero

    specHP = specHPP
    specHP = specHP / 50

    #ppeak2 = getPeaks(specHP, nm, wSize, pkStart, pkStop)   # get peaks after highpass filter


    #plt.plot(nm, Afilter)
    # plt.plot(nm,ppeak2/2)  # show position of peaks 

    aClip = np.copy(A)
    idxStart = nm.searchsorted(pkStart, 'right') - 1  # last index before nm range of interest
    idxStop = nm.searchsorted(pkStop, 'right')  # first index after nm range of interest
    aClip[0:idxStart] = 0
    aClip[idxStop:] = 0

    pkI, _ = find_peaks(aClip, height=5, prominence=4, distance=10)
    # pkI = ppeak2.nonzero()[0]
    pk = As[pkI] + 1.05    # offset in graph Y units, to place above spectrum plot line, not on it
    nmPk = nm[pkI]
    plt.plot(nmPk, pk, "kv")
    for i in pkI:
        x = nm[i]
        y = As[i] + 2.1
        s = ("%5.1f" % x)
        if (x > 420):
            plt.text(x,y,s,horizontalalignment='center')

    printPeaks(pkI, aClip, nm, pkStart, pkStop)

    # ==== end find peaks ======================


    plt.draw()
    plt.pause(0.001)
    
    minAmp = A.min()
    maxAmp = A.max()
    avgAmp = A.mean()
    print("%5.1f, %5.1f, %5.1f" % (minAmp, avgAmp, maxAmp))
    timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    fnameOut = os.path.join(outDir, timeStamp + '_spec.csv')
    df = pd.DataFrame({'nm':nm, 'counts':A, 'cSqrt':As})
    df.to_csv(fnameOut, float_format="%5.2f", sep=',', index=None) # write output data


    time.sleep(0.1)

    #for i in range(0,2048,100):
    #    print("%5.1f, %5.1f" % (nm[i],A[i]))

