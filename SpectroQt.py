# Read and display spectrum
# thanks to:
# https://stackoverflow.com/questions/12459811/how-to-embed-matplotlib-in-pyqt
# https://pythonpyqt.com/qtimer/ 
#
# J.Beale 5/26/2024

import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from PyQt5.QtWidgets import QWidget,QListWidget,QInputDialog,QGridLayout,QLabel
from PyQt5.QtCore import QSize,QTimer,QDateTime

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from seabreeze.spectrometers import Spectrometer
import bottleneck as bn  # boxcar filter
from scipy.signal import find_peaks

import random
import numpy as np
import numpy.polynomial.polynomial as poly # fit to calibration curve


# get array of wavelengths in nm for each sensor pixel
def get_nm():
    # y = Ac + Bx + Cx^2 + Dx^3
    coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
    spIdx = np.linspace(0, specSensorPixels-1, specSensorPixels)
    nm = poly.polyval(spIdx, coefs)
    return (nm)

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
def printPeaks(pkIdx, A, pos, start, stop):
    dCount = 0
    current_time=QDateTime.currentDateTime()
    tstamp = current_time.toString('yy-MM-dd hh:mm:ss')

    for i in pkIdx:  # display position of peaks in nm
        pk = pos[i]
        amp = A[i]
        if (pk >= start) and (pk <= stop):
            print("%s  %d : %5.3f nm  %5.3f counts" % (tstamp, i,pk, amp))
            dCount += 1
    # print("Peaks: %d" % dCount)

class Window(QDialog):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)
        self.setWindowTitle("SpectroQt")
        self.setMinimumSize(QSize(1500, 500))

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.doScan = QPushButton('1-Scan')
        self.doScan.clicked.connect(self.plot)
        self.doBase = QPushButton('Baseline')
        self.doBase.clicked.connect(self.doBaseline)
        self.bRef = QPushButton('Blank-Ref')
        self.bRef.clicked.connect(self.doSetReference)
        self.Reset = QPushButton('Reset')
        self.Reset.clicked.connect(self.doReset)

        self.startBt=QPushButton('Start')
        self.startBt.clicked.connect(self.startTimer)
        self.endBtn=QPushButton('Stop')
        self.endBtn.clicked.connect(self.endTimer)
        self.yScaleA=QPushButton('Y rescale')
        self.yScaleA.clicked.connect(self.yRescale)
        self.setExp=QPushButton('Exposure')
        self.setExp.clicked.connect(self.setExposure)
        self.setPk=QPushButton('Peaks')
        self.setPk.clicked.connect(self.togglePeaks)

        #self.set = QLineEdit()

        self.header=QLabel('Spectrometer App')
        self.status=QLabel('Scan Paused')

        # set the layout
        #layout = QVBoxLayout()
        layout=QGridLayout()

        layout.addWidget(self.header, 0,0)
        layout.addWidget(self.status, 0,1,1,1)
        layout.addWidget(self.toolbar,1,0)
        layout.addWidget(self.canvas, 2,0,5,9) # main plot goes here
        layout.addWidget(self.startBt,7,0,1,1)
        layout.addWidget(self.endBtn, 7,1,1,1)
        layout.addWidget(self.doScan, 7,2,1,1)
        layout.addWidget(self.doBase, 7,3,1,1)
        layout.addWidget(self.bRef,   7,4,1,1)
        layout.addWidget(self.Reset,  7,5,1,1)
        layout.addWidget(self.yScaleA,7,6,1,1)
        layout.addWidget(self.setExp, 7,7,1,1)
        layout.addWidget(self.setPk,  7,8,1,1)

        self.setLayout(layout)

        self.timer=QTimer()
        self.timer.timeout.connect(self.showTime)

        self.exposure_ms = 8.0 # spectrometer integration time in milliseconds        
        self.spec = Spectrometer.from_first_available()    # only one Ocean Optics spectrometer?
        self.spec.integration_time_micros(int(self.exposure_ms * 1000))
        self.averages = 100  # how many scans to average together
        self.boxcar = 5    # total adjacent pixels to boxcar-average (must be odd)        
        self.yRange = (0,specSensorMaxVal)  # vertical scaling of plot
        self.xRange = (365,900)  # horizontal axis range in nm
        self.baseline = np.zeros(specSensorPixels)
        self.blankRef = np.zeros(specSensorPixels)
        self.refSet = False  # haven't taken a blank reference yet
        self.yLabelText = 'intensity (counts)'

        self.showPeaks = False

    # replace all NaN values with 0
    def rmN(self, A):
        A[np.isnan(A)] = 0
        return(A)

    # boxcar-filter array A with window size w
    # w should be an odd number, otherwise result has offset
    def filt(self, A, w):
        b = bn.move_mean(A, window=w) # should left-shift result by (n-1)/2
        b = self.rmN(b)
        res = np.roll(b, -int((w-1)/2) )
        return(res)

    # get a spectrum which is the average of N acquisitions, boxcar Avg pixels together
    def getSpec(self, N, Avg):    
        A = self.spec.intensities()
        for i in range(N-1):
            A += self.spec.intensities()
        A = A / N  # convert the sum to an average
        if (Avg > 1): # average together neighboring pixels?
            A = self.filt(A, Avg)
        A = A - self.baseline  # remove the dark frame            
        return(A)

    def drawPlot(self):
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(nm, self.A)
        self.ax.grid('both')
        self.ax.set_ylim(self.yRange)
        self.ax.set_xlim(self.xRange)
        self.ax.set_ylabel(self.yLabelText)
        self.ax.set_xlabel("nm")
        minor_locator = AutoMinorLocator(5)
        self.ax.xaxis.set_minor_locator(minor_locator)

        self.canvas.draw() # display graph on Qt canvas 

    def plot(self):
        ''' display spectrum data on Qt canvas '''
        self.A = self.getSpec(self.averages, self.boxcar)

        if (self.refSet):
            self.A = np.clip(self.A,0.00001,4095)
            ratio = np.clip((self.blankRef/self.A),0.00001,1000)
            self.A = np.log10(ratio)
        self.drawPlot()
        if (self.showPeaks):
            self.doPeaks()

    # set the sensor baseline (dark frame). All values should be on [0,4095] interval
    def doBaseline(self):
        self.baseline = np.zeros(specSensorPixels)
        self.baseline = self.getSpec(self.averages, self.boxcar)

    # get reference spectrum, switch to displaying absorbance in OD
    def doSetReference(self):
        self.blankRef = np.zeros(specSensorPixels)
        self.blankRef = self.getSpec(self.averages, self.boxcar)
        self.blankRef = np.clip(self.blankRef,0.00001,4095)
        self.refSet = True
        self.yRange = (0, 2.5)  # appropriate units for density
        self.yLabelText = 'Absorbance (OD)'

    # exit absorption mode, return to normal spectrum intensity
    def doReset(self):
        self.refSet = False
        self.yRange = (0,specSensorMaxVal)  # reset vertical scaling of plot
        self.yLabelText = 'intensity (counts)'

    def yRescale(self):               
        yTop = np.max( self.filt(self.A, 15) ) * 1.1
        self.yRange = (0,yTop)  # reset Y axis range        
        self.drawPlot()

    def setExposure(self):
        exp, done1 = QInputDialog.getInt(
            self, 'Exposure Setting', 'New exposure (ms):') 
        if done1:
            exp = np.clip(exp,3,10000)  # USB2000 minimum exp. is 3 msec
            self.exposure_ms = exp
            self.spec.integration_time_micros(int(self.exposure_ms * 1000))            
    def doPeaks(self):
        pkI = getPeaks(self.A, pkStart, pkStop, pkHeight, pkProminence, pkDistance)
        printPeaks(pkI, self.A, nm, pkStart, pkStop)

    def togglePeaks(self):
        self.showPeaks = not self.showPeaks

    def showTime(self):
        current_time=QDateTime.currentDateTime()
        formatted_time=current_time.toString('yyyy-MM-dd hh:mm:ss')
        self.status.setText(formatted_time)
        self.plot()

    def startTimer(self):
        self.timer.start(3000)
        self.startBt.setEnabled(False)
        self.endBtn.setEnabled(True)

    def endTimer(self):
        self.timer.stop()
        self.startBt.setEnabled(True)
        self.endBtn.setEnabled(False)

# ================================================================================

if __name__ == '__main__':
    app = QApplication(sys.argv)

    # set gobal vars related to spectrometer sensor hardware
    specSensorPixels = 2048  # count of sensor pixels
    specSensorMaxVal = 4095  # max possible 12-bit sensor value
    pkStart=375 # find no peaks below this
    pkStop=780  # find no peaks above this
    pkHeight=5  # peak finding parameters
    pkProminence=10
    pkDistance=10

    nm = get_nm()  # global variable with spectrometer calibration array in nanometers
    main = Window()
    main.show()

    sys.exit(app.exec_())
    
