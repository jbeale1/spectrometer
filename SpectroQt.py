# Read and display data from Ocean Optics spectrometer
# thanks to:
# https://stackoverflow.com/questions/12459811/how-to-embed-matplotlib-in-pyqt
# https://pythonpyqt.com/qtimer/ 
# https://physics.nist.gov/PhysRefData/Handbook/Tables/mercurytable2_a.htm
#
# J.Beale 5/27/2024

import sys, os
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QHBoxLayout,QSizePolicy
from PyQt5.QtWidgets import QWidget,QListWidget,QInputDialog,QLabel
from PyQt5.QtCore import QSize,QTimer,QDateTime

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib import rc
from seabreeze.spectrometers import Spectrometer
import bottleneck as bn  # boxcar filter
from scipy.signal import find_peaks

import numpy as np
import numpy.polynomial.polynomial as poly # fit to calibration curve
import pandas as pd

# get array of wavelengths in nm for each sensor pixel
def get_nm():
    # y = Ac + Bx + Cx^2 + Dx^3
    coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
    spIdx = np.linspace(0, specSensorPixels-1, specSensorPixels)
    nm = poly.polyval(spIdx, coefs)
    return (nm)

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

# get data from spectrometer
# average together N acquisitions, and boxcar Avg adjacent pixels together
def getSpec(spec, baseline, N, Avg):    
    A = spec.intensities()
    for i in range(N-1):
        A += spec.intensities()
    A = A / N  # convert the sum to an average
    if (Avg > 1): # average together neighboring pixels?
        A = filt(A, Avg)
    A = A - baseline  # remove the dark frame            
    return(A)

# find index of peaks in spectrum A
def getPeaks(A, pkStart, pkStop, height, prominence, distance):
    aClip = np.copy(A)    
    idxStart = nm.searchsorted(pkStart, 'right') - 1  # last index before nm range of interest
    idxStop = nm.searchsorted(pkStop, 'right')  # first index after nm range of interest
    aClip[0:idxStart] = 0
    aClip[idxStop:] = 0
    pkI, _ = find_peaks(aClip, height=height, 
                        prominence=prominence, distance=distance)
    return (pkI)

# display peak location
def printPeaks(pkIdx, A, pos, start, stop, refSet):
    dCount = 0
    current_time=QDateTime.currentDateTime()
    tstamp = current_time.toString('yy-MM-dd hh:mm:ss')

    for i in pkIdx:  # display position of peaks in nm
        pk = pos[i]
        amp = A[i]
        if refSet:
            units = "OD"
        else:            
            units = "counts"
        if (pk >= start) and (pk <= stop):
            print("%s  %d : %5.3f nm  %5.3f %s" % (tstamp, i,pk, amp, units))
            dCount += 1
    # print("Peaks: %d" % dCount)

# plot signal with labelled peaks ====================
def plotPeaks(fig, ax, nm, A, pkI, timeStamp):

    #fig, ax = plt.subplots()    
    #current_time=QDateTime.currentDateTime()
    #timeStamp=current_time.toString('yyyy-MM-dd hh:mm:ss')

    numFont = {'size': 7}
    rc('font', **numFont)

    yLim = ax.get_ylim()  # actual y limit of axes
    yLim1 = yLim[1] + 0.05*(yLim[1]-yLim[0]) # a little margin
    ax.set_ylim(yLim[0], yLim1)
    yscale = (yLim[1]-yLim[0])/45
    for i in pkI:
        x = nm[i]
        ym = A[i] + yscale*1.05 # offset in Y units to show maker above spectrum plot line, not on it
        yOffset = 0.2
        #if (abs(577-x) < 1.5):
        #    yOffset = 2;
        y1 = A[i] + yscale * ((2.1 * 1.4) + yOffset)
        y2 = A[i] + yscale * ((2.1 * 0.82) + yOffset)
        s1 = ("%5.1f" % x)  # wavelength in nm
        s2 = ("%d" % i)     # sensor pixel; index number
        if (x > pkStart) and (x < pkStop):
            ax.scatter(x, A[i] + yscale*1.05, s=40, marker="v", facecolors='none', 
                       edgecolors='#0000a0') # peak marker
            ax.text(x,y1,s1,horizontalalignment='center') # peak text
            ax.text(x,y2,s2,horizontalalignment='center') # second line of text

    labelFont = {'size': 12}
    rc('font', **labelFont)

    if (pkI.size > 0):
        j = np.argmax(A[pkI])
        pkAmp = A[pkI[j]]   # amplitude of largest peak
        pknm = nm[pkI[j]]   # position of peak in nm
        plt.annotate("%s\n%5.1f nm  (%5.1f)" % (timeStamp,pknm,pkAmp), xy=(0.84,0.88),xycoords='axes fraction')
    else:
        plt.annotate("%s" % (timeStamp), xy=(0.84,0.92),xycoords='axes fraction')                
    # ax.annotate("%s" % (timeStamp), xy=(0.87,0.94),xycoords='axes fraction')
    ax.grid(visible=True, axis='both', color=(0.5, 0.5, 0.5, 0.2),
            linestyle='solid', linewidth='0.5')

# save output files and graph image
def saveData(A,pkI,nm,outDir,labelRaw,timeStamp):
    # timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())

    label = labelRaw.replace(" ", "_")
    label = ''.join(e for e in label if (e.isalnum() or e=="_"))

    fnameOut = os.path.join(outDir, timeStamp + '_spec_' + label + '.csv')
    peakNameOut = os.path.join(outDir, timeStamp + '_peak_' + label + '.csv')    
    plotNameOut = os.path.join(outDir, timeStamp + '_plot_' + label + '.png')

    df = pd.DataFrame({'nm':nm, 'counts':A})
    df.to_csv(fnameOut, float_format="%5.2f", sep=',', index=None) # write waveform data
    dfPk = pd.DataFrame({'index':pkI, 'nm':nm[pkI], 'counts':A[pkI]})
    dfPk.to_csv(peakNameOut, float_format="%5.2f", sep=',', index=None) # write peak data
    plt.savefig(plotNameOut, bbox_inches='tight')

    print("Saved %s" % fnameOut)
    print("Saved %s" % peakNameOut)
    print("Saved %s" % plotNameOut)


# ================================================================================

class Window(QDialog):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("SpectroQt")
        self.setMinimumSize(QSize(1500, 700))
        self.sampleName = ""

        self.figure = plt.figure() # a figure instance to plot on
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.doScan = QPushButton('1-Scan')
        self.doScan.clicked.connect(self.plot1)
        self.doBase = QPushButton('Dark Frame')
        self.doBase.clicked.connect(self.doBaseline)
        self.bRef = QPushButton('Set Ref')
        self.bRef.clicked.connect(self.doSetReference)
        self.Reset = QPushButton('Clear Ref')
        self.Reset.clicked.connect(self.doReset)

        self.startBt=QPushButton('Start')
        self.startBt.clicked.connect(self.startTimer)        
        self.endBtn=QPushButton('Stop')
        self.endBtn.clicked.connect(self.endTimer)
        self.yScaleA=QPushButton('Y rescale')
        self.yScaleA.clicked.connect(self.yRescale)
        self.setExp=QPushButton('Exposure')
        self.setExp.clicked.connect(self.setExposure)
        self.setAvg=QPushButton('Averages')
        self.setAvg.clicked.connect(self.setAverages)
        self.setSm=QPushButton('Smooth')        
        self.setSm.clicked.connect(self.setSmooth)        
        self.setPk=QPushButton('Peaks')
        self.setPk.clicked.connect(self.togglePeaks)
        self.setName=QPushButton('sName')
        self.setName.clicked.connect(self.enterName)
        self.writeBt=QPushButton('save')
        self.writeBt.clicked.connect(self.writeCSV)

        self.header=QLabel('Spectrometer App')
        self.status=QLabel('Scan Paused')
        self.status.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)

        btn1Layout = QHBoxLayout()  # a row of buttons
        btn1Layout.addWidget(self.startBt)
        btn1Layout.addWidget(self.endBtn)
        btn1Layout.addWidget(self.doScan)
        btn1Layout.addWidget(self.doBase)
        btn1Layout.addWidget(self.bRef)
        btn1Layout.addWidget(self.Reset)
        btn1Layout.addWidget(self.yScaleA)
        btn1Layout.addWidget(self.setExp)
        btn1Layout.addWidget(self.setAvg)
        btn1Layout.addWidget(self.setSm)
        btn1Layout.addWidget(self.setPk)
        btn1Layout.addWidget(self.setName)
        btn1Layout.addWidget(self.writeBt)

        hdrLayout = QHBoxLayout()  # header row
        hdrLayout.addWidget(self.toolbar)
        hdrLayout.addWidget(self.status) # current program status      
        # without size policy, QLabel will expand, instead of graph canvas      
        self.status.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        layout = QVBoxLayout()
        layout.addLayout(hdrLayout)
        layout.addWidget(self.canvas)
        layout.addLayout(btn1Layout)
        self.setLayout(layout)

        self.timer=QTimer()
        self.timer.timeout.connect(self.showTime)

        self.exposure_ms = 8.0 # spectrometer integration time in milliseconds        
        self.averages = 20  # how many scans to average together
        self.boxcar = 5    # total adjacent pixels to boxcar-average (must be odd)        
        self.spec = Spectrometer.from_first_available()    # only one Ocean Optics spectrometer?
        self.spec.integration_time_micros(int(self.exposure_ms * 1000))
        self.doReset()       # reset various parameters to intensity, not absorption mode
        self.xRange = (380,1000)  # horizontal axis range in nm
        self.baseline = np.zeros(specSensorPixels)
        self.blankRef = np.zeros(specSensorPixels)
        self.showPeaks = True
        qtime=QDateTime.currentDateTime()
        self.timestamp = qtime.toString('yyyyMMdd_hhmmss')
        self.tString=qtime.toString('yyyy-MM-dd hh:mm:ss')

        self.plot(True)  # draw the first plot

    def drawPlot(self,tstring):
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.figure.subplots_adjust(left=0.07, right=0.97, top=0.9, bottom=0.1)
        self.ax.plot(nm, self.A)
        self.ax.grid('both')
        self.ax.set_ylim(self.yRange)
        self.ax.set_xlim(self.xRange)
        self.ax.set_ylabel(self.yLabelText)
        self.ax.set_xlabel("wavelength  nm")
        self.ax.xaxis.set_major_locator(MultipleLocator(50))
        minor_locator = AutoMinorLocator(5)
        self.ax.xaxis.set_minor_locator(minor_locator)
        self.ax.set_title(self.sampleName)
        if (self.showPeaks):
            plotPeaks(self.figure, self.ax, nm, self.A, self.pkI, tstring)

        self.canvas.draw() # display graph on Qt canvas 

    def plot1(self):
        self.plot(True)
    
    def plot(self, newDat):
        ''' display spectrum data on Qt canvas '''
        
        if (newDat):
            self.status.setText("Running...")        
            self.status.update        
            self.A = getSpec(self.spec, self.baseline, self.averages, self.boxcar) # actually read spectrometer
            
            current_time=QDateTime.currentDateTime() # update timestamp display
            self.tString=current_time.toString('yyyy-MM-dd hh:mm:ss')
            statusString = (" (%d-%d)  t=%dms av=%d sm=%d : %s" % 
                        (self.xRange[0], self.xRange[1], 
                         self.exposure_ms, self.averages, self.boxcar, self.tString))
            self.status.setText(statusString)        

            if (self.refSet):
                self.A = np.clip(self.A,0.00001,4095)
                ratio = np.clip((self.blankRef/self.A),0.00001,1000)
                self.A = np.log10(ratio)
            if (self.showPeaks):
                self.doPeaks()
                
        self.drawPlot(self.tString)

    # set the sensor baseline (dark frame). All values should be on [0,4095] interval
    def doBaseline(self):
        baseZ = np.zeros(specSensorPixels)
        self.baseline = getSpec(self.spec, baseZ, self.averages, self.boxcar) # actually read spectrometer        
        self.plot(True)  # do another acquisition to refresh plot with new baseline


    # get reference spectrum, switch to displaying absorbance in OD
    def doSetReference(self):        
        self.blankRef = getSpec(self.spec, self.baseline, self.averages, self.boxcar) # actually read spectrometer        
        self.blankRef = np.clip(self.blankRef,0.00001,4095)
        self.refSet = True
        self.yRange = (0, 2.5)  # appropriate units for density
        self.yLabelText = 'Absorbance (OD)'
        self.pkHeight=.01
        self.pkProminence=0.001
        self.pkDistance=40

    # exit absorption mode, return to normal spectrum intensity
    def doReset(self):
        self.refSet = False
        self.yRange = (0,specSensorMaxVal)  # reset vertical scaling of plot
        self.yLabelText = 'intensity (counts)'
        self.pkHeight=5  # peak finding parameters in amplitude mode
        self.pkProminence=10
        self.pkDistance=10

    def yRescale(self):     
        Amod = self.A[30:] # first several elements are too noisy
        if (self.refSet):
            yTop = np.max( filt(Amod, 15) ) * 1.1
        else:            
            yTop = np.max( Amod ) * 1.1
        self.yRange = (0,yTop)  # reset Y axis range        
        self.drawPlot()

    def setExposure(self):
        exp, done1 = QInputDialog.getInt(
            self, 'Exposure Setting', 'New exposure (ms):', int(self.exposure_ms)) 
        if done1:
            exp = np.clip(exp,3,10000)  # USB2000 minimum exp. is 3 msec
            self.exposure_ms = exp
            self.spec.integration_time_micros(int(self.exposure_ms * 1000))   

    def setAverages(self):
        avg, done = QInputDialog.getInt(
            self, 'Set Averages', 'Readings to average:', self.averages) 
        if done:
            avg = np.clip(avg,1,10000)  # don't get too crazy
            self.averages = avg            

    def setSmooth(self):
        avg, done = QInputDialog.getInt(
            self, 'Set Smoothing', 'Pixels to boxcar-average:', self.boxcar) 
        if done:
            avg = np.clip(avg,1,252)  # don't get too crazy
            if (avg % 2)==0:  # should not be an even n umber
                avg = avg-1
            self.boxcar = avg            

    def enterName(self):
        self.sampleName, done1 = QInputDialog.getText(
            self, 'Enter Sample Name', 'Sample name:') 
        self.plot(False)  # do another acquisition to refresh plot

    def doPeaks(self):
        self.pkI = getPeaks(self.A, pkStart, pkStop,
                       self.pkHeight, self.pkProminence, self.pkDistance)
        printPeaks(self.pkI, self.A, nm, pkStart, pkStop, self.refSet)

    def togglePeaks(self):
        self.showPeaks = not self.showPeaks

    def showTime(self):
        current_time=QDateTime.currentDateTime()
        formatted_time=current_time.toString('yyyy-MM-dd hh:mm:ss')
        self.timestamp = current_time.toString('yyyyMMdd_hhmmss')
        self.status.setText(formatted_time)
        self.plot1()

    def startTimer(self):
        # calculate how much time to acquire and display one frame
        delay = 150 + int((20 + self.exposure_ms) * self.averages)
        self.timer.start(delay)
        self.startBt.setEnabled(False)
        self.endBtn.setEnabled(True)

    def endTimer(self):
        self.timer.stop()
        self.startBt.setEnabled(True)
        self.endBtn.setEnabled(False)

    def writeCSV(self):
        print("Write csv")
        # timeStamp = time.strftime("%Y%m%d_%H%M%S", time.localtime())
        saveData(self.A,self.pkI,nm,outDir,self.sampleName,self.timestamp)

# ================================================================================

if __name__ == '__main__':
    app = QApplication(sys.argv)

    # set gobal vars related to spectrometer sensor hardware
    specSensorPixels = 2048  # count of sensor pixels
    specSensorMaxVal = 4095  # max possible 12-bit sensor value
    pkStart=380 # find no peaks below this
    pkStop=800  # find no peaks above this
    pkHeight=5  # peak finding parameters in amplitude mode
    pkProminence=10
    pkDistance=20
    outDir = r"C:\Users\beale\Documents\OceanSpec\RawData"

    nm = get_nm()  # global variable with spectrometer calibration array in nanometers
    main = Window()
    main.show()

    sys.exit(app.exec_()) # forwards any PyQt exit error code to calling process
    
