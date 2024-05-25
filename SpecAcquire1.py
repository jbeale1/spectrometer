# acquire spectrum from Ocean Optics USB200
# and find peaks
# 24-May-2024 J.Beale

# https://www.researchgate.net/figure/Fluorescence-spectra-of-chlorophyll-a-l-max-650-nm-and-chlorophyll-b-l-max-670-nm_fig1_272366292

import numpy as np
from matplotlib import pyplot as plt 
from matplotlib import rc
# from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.ticker import AutoMinorLocator
from matplotlib import cm

from seabreeze.spectrometers import Spectrometer
import time  # for time-date and real-time delay
import bottleneck as bn  # boxcar filter
import os  # combine dir + filename into path
import pandas as pd  # write .csv output file
from scipy.signal import find_peaks

import numpy.polynomial.polynomial as poly # fit to calibration curve

outDir = "/home/john/Documents/Spectrometer"  # directory for output data files
titleString = "USB2000 : 6-inch F4T5-CW fluorescent bulb" # plot title

integrationTime = 15000 # spectrometer integration time in microseconds
averages = 400          # how many spectra to average together
cycles = 5            # total number of datasets to produce

boxcar = 1             # boxcar-filter this many adjacent pixels (always odd)
# xRange = (360, 1000)    # display this range in nm

#yRange = (0, 4200)      # this range in intensity units (12-bit ADC)
yRange = (0, 50)      # plot this range in intensity units (12-bit ADC)
pltWidth = 22       # plot output window size, in inches?
pltHeight = 8

#  peak-finding parameters
wSize= 11  # rolling average window size
pkSearchWin = 41 # array element range, must be odd
pkStart = 359  # minimum wavelength peak of interest
pkStop = 1000   # maximum wavelength peak of interest
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
            print("%d : %5.3f nm  %5.3f counts" % (i,pk, amp))
            dCount += 1
    print("Peaks: %d" % dCount)
    
# ================================================================================    
# literature values for neon emission (in nm) that a neon bulb can make    
# "Wavelengths, energy level classifications, and energy levels for the 
# spectrum of neutral neon" Jan.2004 Journal of Physical and Chemical Reference Data
neon_nm = np.array([540.05618, 585.24879, 614.30626, 640.2248,
    667.82762, 703.24131, 724.51666, 837.7608, 878.06226    ])

old_nm = np.array([538.87, 583.99, 613.03, 638.86, 666.37, 701.87,
                   723.06, 836.29, 876.49 ])

testIdx = np.array([533,668,756,835,920,1031,1098,1467,1603 ])

# Literature values for selected argon emission lines (in nm)
# https://physics.nist.gov/PhysRefData/Handbook/Tables/argontable2.htm
argon_nm = np.array([696.5431, 738.3980, 750.3869, 763.5106, 811.5311, 842.4648, 912.2967 ])

# selected mercury emission lines, in nm
# https://www.rp-photonics.com/standard_spectral_lines.html
# https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2_a.htm
mercury_nm  = np.array([365.0153, 404.6563, 435.8328,  546.0735, 576.9598, 579.0663 ])
mercury_idx = np.array([26, 138,  228, 550,  643, 649])

cal_nm = np.array([365.0153, 404.6563, 435.8328, 540.05618, 546.0735, 585.24879, 614.30626, 640.2248,
    667.82762, 703.24131, 724.51666, 837.7608, 878.06226    ])

cal_idx = np.array([26, 138, 228, 533, 550, 668,756,835,920,1031,1098,1467,1603 ])

def doPfit(x,y,n):
    coefs = poly.polyfit(x, y, n)
    x_new = np.linspace(x[0], x[-1], num=len(x)*500)
    ffit = poly.polyval(x_new, coefs)
    plt.plot(x_new, ffit, ".")
    plt.plot(x,y,"+")
    eMax = 0 # maximum error
    for i in range(len(x)):
        xp = x[i] # original sensor pixel index
        yp = y[i] # wavelength in nm
        idx = (np.abs(x_new - xp)).argmin() # index of closest in ffit
        error = ffit[idx]-yp
        if abs(error)>abs(eMax):
            eMax = error
        print("%d, %d, %5.3f, %5.3f, %3.3f" % (xp, idx, yp, ffit[idx], error))
    print("Max error magnitude: %5.3f" % eMax)
    print("coefs: ", coefs)
    plt.show()

# get array of wavelengths in nm for each sensor pixel
def get_nm():
    # y = Ac + Bx + Cx^2 + Dx^3
    #coefs = np.array([3.54859155e+02, 3.56950640e-01, -1.72680824e-05, -1.12499513e-09])
    coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
    spIdx = np.linspace(0, 2047, 2048)
    nm = poly.polyval(spIdx, coefs)
    return (nm)

# ==================================================================
# main code starts here

#doPfit(testIdx, neon_nm,3)
#doPfit(testIdx, old_nm,3)
#doPfit(cal_idx, cal_nm, 3)
#exit()

spec = Spectrometer.from_first_available()
spec.integration_time_micros(integrationTime)

#nm = spec.wavelengths()
nm = get_nm()

# plt.axis([365,1000,0,4096])
#plt.ion()
#plt.show()

print("\nStarting baseline... ",end='')
baseline = getSpec(int(1+averages/2), boxcar)  # get baseline reference
print(" done. Running acquisition...")

for j in range(cycles):

    A = getSpec(averages, boxcar) - baseline
    As = np.sqrt(np.clip(A,0,4095))
    #plt.clf()  # clear current plot figure
    
    if plt.get_fignums(): # empty list is false
        plt.close()

    fig, ax = plt.subplots()

    timeStamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    # fig, ax = plt.subplots()
    # plt.figure(figsize=(pltWidth, pltHeight))
    fig.set_figheight(pltHeight)
    fig.set_figwidth(pltWidth)
    

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)

    ax.set_xlim(xRange)
    ax.set_ylim(yRange)
    ax.set_title(titleString)
    ax.set_xlabel('wavelength, nm')
    ax.set_ylabel('sqrt(intensity)')
    ax.plot(nm, As, linewidth=0.6, color='#000080') # show the spectrum plot


    # find peaks ============================


    aClip = np.copy(A)
    idxStart = nm.searchsorted(pkStart, 'right') - 1  # last index before nm range of interest
    idxStop = nm.searchsorted(pkStop, 'right')  # first index after nm range of interest
    aClip[0:idxStart] = 0
    aClip[idxStop:] = 0

    # pkI, _ = find_peaks(aClip, height=2, prominence=8, distance=50)
    pkI, _ = find_peaks(aClip, height=5, prominence=10, distance=3)
    # pkI = ppeak2.nonzero()[0]
    #pk = As[pkI] + 1.05    # offset in graph Y units, to place above spectrum plot line, not on it
    #nmPk = nm[pkI]
    # plt.plot(nmPk, pk, "kv")
    numFont = {'size': 7}
    rc('font', **numFont)

    for i in pkI:
        x = nm[i]
        ym = As[i] + 1.05 # offset in Y units to show maker above spectrum plot line, not on it
        yOffset = 0.2
        #if (abs(390-i) < 3):
        #    continue
        if (abs(577-x) < 1.5):
            yOffset = 2;
        y1 = As[i] + (2.1 * 1.2) + yOffset
        y2 = As[i] + (2.1 * 0.85) + yOffset
        s1 = ("%5.1f" % x)
        s2 = ("nm")
        # s1 = ("%d" % i)
        # s1 = ("%d" % i)
        if (x > pkStart):
            ax.scatter(x, As[i] + 1.05, s=40, marker="v", facecolors='none', 
                       edgecolors='#0000a0') # peak marker
            ax.text(x,y1,s1,horizontalalignment='center') # peak text
            ax.text(x,y2,s2,horizontalalignment='center')

    printPeaks(pkI, aClip, nm, pkStart, pkStop)

    # ==== end find peaks ======================
    labelFont = {'size': 12}
    rc('font', **labelFont)

    # plt.annotate("%s\n%d peaks" % (timeStamp,pkI.size), xy=(0.84,0.92),xycoords='axes fraction')
    ax.annotate("%s" % (timeStamp), xy=(0.87,0.94),xycoords='axes fraction')
    ax.grid(visible=True, axis='both', color=(0.5, 0.5, 0.5, 0.2),
            linestyle='solid', linewidth='0.5')


    timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    fnameOut = os.path.join(outDir, timeStamp + '_spec.csv')
    peakNameOut = os.path.join(outDir, timeStamp + '_peak.csv')
    plotNameOut = os.path.join(outDir, timeStamp + '_plot.png')

    #plt.axes().xaxis.set_major_locator(MultipleLocator(50))
    #plt.axes().xaxis.set_minor_locator(MultipleLocator(5))

    #fig.tight_layout() # call after all axes drawn
    # plt.tight_layout()
    plt.savefig(plotNameOut, bbox_inches='tight')
    plt.draw()
    plt.pause(0.001)
    
    df = pd.DataFrame({'nm':nm, 'counts':A, 'cSqrt':As})
    df.to_csv(fnameOut, float_format="%5.2f", sep=',', index=None) # write waveform data
    dfPk = pd.DataFrame({'index':pkI, 'nm':nm[pkI], 'counts':A[pkI]})
    dfPk.to_csv(peakNameOut, float_format="%5.2f", sep=',', index=None) # write peak data

    minAmp = A.min()
    maxAmp = A.max()
    avgAmp = A.mean()
    print("%s , %5.1f, %5.1f, %5.1f" % (fnameOut, minAmp, avgAmp, maxAmp))

    time.sleep(0.1)
