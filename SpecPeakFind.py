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

outDir = "/home/john/Documents/Spectrometer"  # directory for output data files

integrationTime = 200000 # spectrometer integration time in microseconds
averages = 50          # how many spectra to average together
boxcar = 5             # boxcar-filter this many adjacent pixels (always odd)
xRange = (360, 1000)    # display this range in nm
xRange = (600, 800)    # display this range in nm
yRange = (0, 4200)      # this range in intensity units (12-bit ADC)

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
    plt.clf()  # clear current plot figure
    plt.plot(nm, A)
    plt.xlim(xRange)
    plt.ylim(yRange)
    plt.title("Spectrum")
    plt.xlabel('wavelength, nm')
    plt.ylabel('intensity')
    plt.grid(axis='both')

    plt.draw()
    plt.pause(0.001)
    
    minAmp = A.min()
    maxAmp = A.max()
    avgAmp = A.mean()
    print("%5.1f, %5.1f, %5.1f" % (minAmp, avgAmp, maxAmp))
    timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    fnameOut = os.path.join(outDir, timeStamp + '_spec.csv')
    
    df = pd.DataFrame({'nm':nm, 'counts':A})
    df.to_csv(fnameOut, float_format="%5.2f", sep=',', index=None) # write output data
    time.sleep(0.1)

    #for i in range(0,2048,100):
    #    print("%5.1f, %5.1f" % (nm[i],A[i]))

