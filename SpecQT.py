import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout

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

        # A button connected to `plot` method
        self.button = QPushButton('Plot')
        self.button.clicked.connect(self.plot)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)

        integrationTime = 100000 # spectrometer integration time in microseconds
        self.spec = Spectrometer.from_first_available()    # only one Ocean Optics spectrometer?
        self.spec.integration_time_micros(integrationTime)
        self.nm = self.get_nm()

    # get array of wavelengths in nm for each sensor pixel
    def get_nm(self):
        # y = Ac + Bx + Cx^2 + Dx^3
        coefs = np.array([ 3.55900745e+02, 3.54281751e-01, -1.50846947e-05, -1.70072117e-09 ])
        spIdx = np.linspace(0, 2047, 2048)
        nm = poly.polyval(spIdx, coefs)
        return (nm)

    def plot(self):
        ''' plot some random stuff '''

        #A = getSpec(averages, boxcar)
        self.A = self.spec.intensities()

        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.plot(self.nm, self.A)
        ax.grid('both')
        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())
    
