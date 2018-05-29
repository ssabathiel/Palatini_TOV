import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QCheckBox, QHBoxLayout
from PyQt5.QtCore import Qt
#import PyQt4
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"


#from PyQt5.QtCore import *
#from PyQt5.QtGui import *

import random
import glob


sil_fR=1
sil_fRQ=1
sil_GR=1
diff=0
profile = 0


###############
## Get data from files to arrays
###############


fR_Rp = glob.glob('Results/TOV_output_fR_*')[0][25:32] #35]
fR_Rp +=   '{' + glob.glob('Results/TOV_output_fR_*')[0][32:35] + '}'
print("Rp=", fR_Rp)

fRQ_Rp = glob.glob('Results/TOV_output_fRQ_*')[0][26:33] #35]
fRQ_Rp +=   '{' + glob.glob('Results/TOV_output_fRQ_*')[0][33:36] + '}'
print("fRQ,Rp=", fRQ_Rp)

fRQ_Rq = glob.glob('Results/TOV_output_fRQ_*')[0][39:46]
fRQ_Rq +=   '{' + glob.glob('Results/TOV_output_fRQ_*')[0][46:49] + '}'
print("fRQ,Rq=", fRQ_Rq)


#fRQ_Rq = glob.glob('Results/TOV_output_fRQ_*')[0][39:49]



with open(glob.glob('Results/TOV_output_fR_*')[0]) as f:
#with open('Results/TOV_output_fR_Rp_2.1*10^-10') as f:
    lines = f.readlines()
    x = [line.split()[2] for line in lines]
    y = [line.split()[1] for line in lines]
    z = [line.split()[0] for line in lines]

with open(glob.glob('Results/TOV_output_fRQ_*')[0]) as f:
    lines = f.readlines()
    x2 = [line.split()[2] for line in lines]
    y2 = [line.split()[1] for line in lines]
    z2 = [line.split()[0] for line in lines]


with open('Results/TOV_output_GR') as f:
    lines = f.readlines()
    x5 = [line.split()[2] for line in lines]
    y5 = [line.split()[1] for line in lines]
    z5 = [line.split()[0] for line in lines]







for i in range(len(x)):
    x[i]=np.log10(float(x[i])/100000.)
    z[i] = np.log10(float(z[i]))
for i in range(len(x5)):
    x5[i]=np.log10(float(x5[i])/100000.)
    z5[i]=np.log10(float(z5[i]))
for i in range(len(x2)):
    x2[i]=np.log10(float(x2[i])/100000.)
    z2[i]=np.log10(float(z2[i]))







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

        # Just some button connected to `plot` method
        self.button = QPushButton('Plot')
        #self.button.clicked.connect(self.plot)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        layout.addWidget(self.button)
        self.setLayout(layout)

        #layout = QHBoxLayout()
        self.btn_GR = QCheckBox("GR")
        self.btn_fR =  QCheckBox("fR")
        self.btn_fRQ = QCheckBox("fRQ")

        self.btn_GR.setChecked(False)
        self.btn_fR.setChecked(False)
        self.btn_fRQ.setChecked(False)

        #self.b1.stateChanged.connect(lambda: self.btnstate(self.b1))
        self.btn_GR.stateChanged.connect(self.plot)
        self.btn_fR.stateChanged.connect(self.plot)
        self.btn_fRQ.stateChanged.connect(self.plot)

        layout.addWidget(self.btn_GR)
        layout.addWidget(self.btn_fR)
        layout.addWidget(self.btn_fRQ)


       # self.plot(self.b1)




        #self.b1.setChecked.connect(self.plot)

    def plot(self):
        fR_label = 'f(R)\n  $(R_p=$' + r'$%s$' % fR_Rp + r'$)$'
        fRQ_label = 'f(R,Q) \n  $(R_p=$' + r'$%s$' % fRQ_Rp + '$)$\n  ' + r'$(R_q=$' + r'$%s$' % fRQ_Rq + '$)$'

        # Plot R-M
        #if sil_fR == 1:
        if self.btn_fR.isChecked() == True:
            plt.plot(x, y, label=fR_label, color='blue')
            plt.title('fR/GR comparison for Mass and Radius for $R_p=2.1*10^{80}$')
            plt.legend()

        if sil_fRQ == 1:
            plt.plot(x2, y2, label=fRQ_label, color='orange')
        if sil_GR == 1:
            plt.plot(x5, y5, label='GR', color='black')
            # TEXT
            plt.xlabel('log R [km]')
            plt.ylabel('$M/M_0$')

            #plt.show()


        ''' plot some random stuff '''
        # random data
        #data = [random.random() for i in range(10)]

        # instead of ax.hold(False)
        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        # discards the old graph
        # ax.hold(False) # deprecated, see above

        # plot data
        if self.btn_GR.isChecked() == True:
            ax.plot(x5, y5, label='GR', color='black')
            ax.legend()

        if self.btn_fR.isChecked() == True:
            ax.plot(x, y, label=fR_label, color='red')
            ax.legend()

        if self.btn_fRQ.isChecked() == True:
            ax.plot(x2, y2, label=fRQ_label, color='blue')
            ax.legend()

        ax.set_title('f(R),f(R,Q),GR comparison for Mass and Radius')
        ax.set_xlim(0.5,4.0)

        # refresh canvas
        self.canvas.draw()




    def btnstate(self, b):
        if b.text() == "fR":
            if b.isChecked() == True:
                self.plot(b)



if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Window()
    main.show()

    sys.exit(app.exec_())