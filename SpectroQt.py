# Read and display spectrum
# thanks to:
# https://stackoverflow.com/questions/12459811/how-to-embed-matplotlib-in-pyqt
# https://pythonpyqt.com/qtimer/ 
#
# J.Beale 5/26/2024

import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout
from PyQt5.QtWidgets import QWidget,QListWidget,QGridLayout,QLabel
from PyQt5.QtCore import QTimer,QDateTime

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from seabreeze.spectrometers import Spectrometer

import random
import numpy as np
import numpy.polynomial.polynomial as poly # fit to calibration curve


class Window(QDialog):

    def __init__(self, parent=None):
        super(Window, self).__init__(parent)

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        self.doScan = QPushButton('Scan')
        self.doScan.clicked.connect(self.plot)

        self.startBtn=QPushButton('Start')
        self.startBtn.clicked.connect(self.startTimer)
        self.endBtn=QPushButton('Stop')
        self.endBtn.clicked.connect(self.endTimer)
        self.yScaleAuto=QPushButton('Y rescale')
        self.yScaleAuto.clicked.connect(self.yRescale)

        self.header=QLabel('Spectrometer App')
        self.status=QLabel('Scan Paused')

        # set the layout
        #layout = QVBoxLayout()
        layout=QGridLayout()

        layout.addWidget(self.header,0,0)
        layout.addWidget(self.status,0,1,1,1)
        layout.addWidget(self.toolbar,1,0)
        layout.addWidget(self.canvas,2,0,5,5)
        layout.addWidget(self.startBtn,7,0,1,1)
        layout.addWidget(self.endBtn,7,1,1,1)
        layout.addWidget(self.doScan,7,2,1,1)
        layout.addWidget(self.yScaleAuto,7,3,1,1)

        self.setLayout(layout)

        self.timer=QTimer()
        self.timer.timeout.connect(self.showTime)


        integrationTime = 100000 # spectrometer integration time in microseconds
        self.spec = Spectrometer.from_first_available()    # only one Ocean Optics spectrometer?
        self.spec.integration_time_micros(integrationTime)
        self.nm = self.get_nm()
        self.yRange = (0,4095)  # vertical scaling of plot

    # get array of wavelengths in nm for each sensor pixel
    def get_nm(self):
        # y = Ac + Bx + Cx^2 + Dx^3
        coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
        spIdx = np.linspace(0, 2047, 2048)
        nm = poly.polyval(spIdx, coefs)
        return (nm)

    def plot(self):
        ''' display spectrum data on Qt canvas '''

        #A = getSpec(averages, boxcar)
        self.A = self.spec.intensities()

        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.plot(self.nm, self.A)
        self.ax.grid('both')
        self.ax.set_ylim(self.yRange)

        self.canvas.draw() # display graph on Qt canvas 

    def yRescale(self):
        yLim = self.ax.get_ylim()  # actual y limit of axes
        yTop = int(np.max(self.A) * 1.1)
        self.yRange = (0,yTop)  # reset Y axis range        

    def showTime(self):
        current_time=QDateTime.currentDateTime()
        formatted_time=current_time.toString('yyyy-MM-dd hh:mm:ss')
        self.status.setText(formatted_time)
        self.plot()

    def startTimer(self):
        self.timer.start(500)
        self.startBtn.setEnabled(False)
        self.endBtn.setEnabled(True)

    def endTimer(self):
        self.timer.stop()
        self.startBtn.setEnabled(True)
        self.endBtn.setEnabled(False)

# ================================================================================

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())
    
