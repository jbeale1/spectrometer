# acquire spectrum from Ocean Optics USB200
# calculate sample absorption and find peaks
# 25-May-2024 J.Beale

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

# outDir = "/home/john/Documents/Spectrometer"  # directory for output data files
outDir = r"C:\Users\beale\Documents\OceanSpec\RawData"

titleString = "Spectrum" # plot title

#integrationTime = 15000 # spectrometer integration time in microseconds
integrationTime = 8000 # spectrometer integration time in microseconds
averages = 500          # how many spectra to average together
cycles = 4            # total number of datasets to produce

boxcar = 5             # boxcar-filter this many adjacent pixels (always odd)
# xRange = (360, 1000)    # display this range in nm

#yRange = (0, 4200)      # this range in intensity units (12-bit ADC)
yRange = (0, 50)      # plot this range in intensity units (12-bit ADC)
pltWidth = 22       # plot output window size, in inches?
pltHeight = 8

#  peak-finding parameters
wSize= 11  # rolling average window size
pkSearchWin = 41 # array element range, must be odd
pkStart = 365  # min wavelength on graph
pkStop = 950   # max wavelength on graph
pkEnd = 830    # don't display any peaks beyond this (nm)
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
    

# get array of wavelengths in nm for each sensor pixel
def get_nm():
    # y = Ac + Bx + Cx^2 + Dx^3
    #coefs = np.array([3.54859155e+02, 3.56950640e-01, -1.72680824e-05, -1.12499513e-09])
    coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
    spIdx = np.linspace(0, 2047, 2048)
    nm = poly.polyval(spIdx, coefs)
    return (nm)

# find index of peaks in spectrum A
def getPeaks(A, height, prominence, distance):
    aClip = np.copy(A)
    idxStart = nm.searchsorted(pkStart, 'right') - 1  # last index before nm range of interest
    idxStop = nm.searchsorted(pkStop, 'right')  # first index after nm range of interest
    aClip[0:idxStart] = 0
    aClip[idxStop:] = 0
    pkI, _ = find_peaks(aClip, height=height, 
                        prominence=prominence, distance=distance)
    return (pkI)

# plot signal with labelled peaks ====================
def plotPeaks(nm, A, pkI, doSqrt, sampleName):
    if plt.get_fignums(): # empty list is false
        plt.close()

    if (doSqrt):
        As = np.sqrt(np.clip(A,0,4095))
        yLabelText = 'sqrt(intensity)'
    else:
        As = A        
        yLabelText = 'Absorbance (OD)'

    fig, ax = plt.subplots()
    timeStamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    fig.set_figheight(pltHeight)
    fig.set_figwidth(pltWidth)

    minor_locator = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.set_xlim(xRange)
    ax.set_title(sampleName)
    ax.set_xlabel('wavelength, nm')
    ax.set_ylabel(yLabelText)
    ax.plot(nm, As, linewidth=0.6, color='#000080') # show the spectrum plot

    numFont = {'size': 7}
    rc('font', **numFont)

    yLim = ax.get_ylim()  # actual y limit of axes
    yLim1 = yLim[1] + 0.05*(yLim[1]-yLim[0]) # a little margin
    ax.set_ylim(yLim[0], yLim1)
    yscale = (yLim[1]-yLim[0])/45
    for i in pkI:
        x = nm[i]
        ym = As[i] + yscale*1.05 # offset in Y units to show maker above spectrum plot line, not on it
        yOffset = 0.2
        #if (abs(577-x) < 1.5):
        #    yOffset = 2;
        y1 = As[i] + yscale * ((2.1 * 1.2) + yOffset)
        y2 = As[i] + yscale * ((2.1 * 0.82) + yOffset)
        s1 = ("%5.1f" % x)  # wavelength in nm
        s2 = ("%d" % i)     # sensor pixel; index number
        if (x > pkStart) and (x < pkEnd):
            ax.scatter(x, As[i] + yscale*1.05, s=40, marker="v", facecolors='none', 
                       edgecolors='#0000a0') # peak marker
            ax.text(x,y1,s1,horizontalalignment='center') # peak text
            ax.text(x,y2,s2,horizontalalignment='center') # second line of text

    labelFont = {'size': 12}
    rc('font', **labelFont)

    # plt.annotate("%s\n%d peaks" % (timeStamp,pkI.size), xy=(0.84,0.92),xycoords='axes fraction')
    ax.annotate("%s" % (timeStamp), xy=(0.87,0.94),xycoords='axes fraction')
    ax.grid(visible=True, axis='both', color=(0.5, 0.5, 0.5, 0.2),
            linestyle='solid', linewidth='0.5')

    plt.draw()
    plt.pause(0.5)


# display basic stats on console
def showStats(A, label):
    minAmp = A.min()
    maxAmp = A.max()
    avgAmp = A.mean()
    print("%s , %5.1f, %5.1f, %5.1f" % (label, minAmp, avgAmp, maxAmp))

# save output files and graph image
def saveData(A,outDir,label):
    As = np.sqrt(np.clip(A,0,4095))

    timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
    fnameOut = os.path.join(outDir, timeStamp + '_spec_' + label + '.csv')
    peakNameOut = os.path.join(outDir, timeStamp + '_peak_' + label + '.csv')
    plotNameOut = os.path.join(outDir, timeStamp + '_plot_' + label + '.png')

    df = pd.DataFrame({'nm':nm, 'counts':A, 'cSqrt':As})
    df.to_csv(fnameOut, float_format="%5.2f", sep=',', index=None) # write waveform data
    dfPk = pd.DataFrame({'index':pkI, 'nm':nm[pkI], 'counts':A[pkI]})
    dfPk.to_csv(peakNameOut, float_format="%5.2f", sep=',', index=None) # write peak data

    plt.savefig(plotNameOut, bbox_inches='tight')
    showStats(A, fnameOut)
    plt.pause(0.5)

# =========================================================================================
#    main code starts here
# =========================================================================================

spec = Spectrometer.from_first_available()    # assume there's only one Ocean Optics spectrometer
spec.integration_time_micros(integrationTime)

# nm = spec.wavelengths()  # built-in saved calibration, not necessarily the best accuracy
nm = get_nm()   # manual calibration

print("\nTaking dark frame... ")
baseline = getSpec(int(1+averages), boxcar)  # get baseline reference
input("Insert empty reference and press enter...")

A = getSpec(averages, boxcar) - baseline    
pkI = getPeaks(A, height=5, prominence=10, distance=3)
printPeaks(pkI, A, nm, pkStart, pkStop)
plotPeaks(nm, A, pkI, True, 'input spectrum')
saveData(A,outDir,'spec')

sampleName = input("Type sample name and press enter...")
if (sampleName == ''):
    sampleName = "Sample Absorbance"

B = getSpec(averages, boxcar) - baseline    
pkI = getPeaks(B, height=5, prominence=10, distance=3)
plotPeaks(nm, B, pkI, True, "sample spectrum")
input("Press enter to see absorption...")

RatioP = np.clip(A/B, 0.00001, 100)
Abs = np.log10(RatioP)
Abs = np.clip(Abs, -0.05, 10) # don't really expect much below 0
pkI = getPeaks(Abs, height=.01, prominence=0.001, distance=40)
printPeaks(pkI, Abs, nm, pkStart, pkStop)
plotPeaks(nm, Abs, pkI, False, sampleName)
saveData(Abs,outDir, sampleName)

input("To exit, oddly enough, press enter...")

# time.sleep(2)
