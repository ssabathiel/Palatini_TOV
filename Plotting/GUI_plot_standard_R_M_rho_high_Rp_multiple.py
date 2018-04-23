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
import subprocess
#subprocess.call("command1")
#subprocess.call(["command1", "arg1", "arg2"])

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

GR_rad = np.empty(0)    #R
fR_rad = np.empty(0)
fRQ_rad = np.empty(0)

GR_m = np.empty(0)
fR_m = np.empty(0)
fRQ_m = np.empty(0)

GR_rho = np.empty(0)
fR_rho = np.empty(0)
fRQ_rho = np.empty(0)



for i in range(len(glob.glob('Results/TOV_output_fR_*'))):
    with open(glob.glob('Results/TOV_output_fR_*')[0]) as f:
        lines = f.readlines()
        np.concatenate((fR_rad, np.asarray( [line.split()[2] for line in lines]  )   ), axis=0 )
        np.concatenate((fR_m,[line.split()[1] for line in lines] )  , axis=0)
        np.concatenate((fR_rho, [line.split()[0] for line in lines] ), axis=0)

for i in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
    with open(glob.glob('Results/TOV_output_fRQ_*')[0]) as f:
        lines = f.readlines()
        np.concatenate((fRQ_rad, [line.split()[2] for line in lines]), axis=0)
        np.concatenate((fRQ_m, [line.split()[1] for line in lines]), axis=0)
        np.concatenate((fRQ_rho, [line.split()[0] for line in lines]), axis=0)


with open('Results/TOV_output_GR') as f:
    lines = f.readlines()
    GR_rad = [line.split()[2] for line in lines]
    GR_m = [line.split()[1] for line in lines]
    GR_rho = [line.split()[0] for line in lines]


print(type(fR_rad))
print(fR_rad.shape)



for j in range(len(glob.glob('Results/TOV_output_fR_*'))):
    for i in range(len(fR_rad[j][:])):
        fR_rad[j][i]=np.log10(float(fR_rad[j][i])/100000.)
        fR_rho[j][i] = np.log10(float(fR_rho[j][i]))

for j in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
    for i in range(len(fRQ_rad[j])):
        fRQ_rad[j][i] = np.log10(float(fRQ_rad[j][i]) / 100000.)
        fRQ_rho[j][i] = np.log10(float(fRQ_rho[j][i]))

for i in range(len(GR_rad)):
    GR_rad[i]=np.log10(float(GR_rad[i])/100000.)
    GR_rad[i]=np.log10(float(GR_rad[i]))


#Cool python tutorials on classes
#https://jeffknupp.com/blog/2014/06/18/improve-your-python-python-classes-and-object-oriented-programming/

#Class Model/Animal inheritates "feature" to subclass fR_Model/Dog or fRQ_Model/Cat
#https://codereview.stackexchange.com/questions/145740/structure-of-inheritance-between-animal-classes-and-methods

#See inheritance with class MyWindow, which can use MyWindow.show() without having it defined: here
#https://stackoverflow.com/questions/3972158/how-to-plot-on-my-gui?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa

#Run .exe from python
#https://www.cyberciti.biz/faq/python-execute-unix-linux-command-examples/

class fR_Model:
    def __init__(self,ID):
        self.ID = ID
        self.rad = np.empty(0)
        self.mass = np.empty(0)
        self.rho = np.empty(0)
        self.Visible = False

    def setData(self, r, m, rh):
        self.rad = r
        self.mass = m
        self.rho = rh


class fRQ_Model:
    def __init__(self, ID):
        self.ID = ID
        self.rad = np.empty(0)
        self.mass = np.empty(0)
        self.rho = np.empty(0)
        self.Visible = False

    def setData(self, r, m, rh):
        self.rad = r
        self.mass = m
        self.rho = rh

    #def setVisible(self, booly):
        #self.Visible = booly


fR_models = np.empty(0)

for j in range(len(glob.glob('Results/TOV_output_fR_*'))):
    np.append(fR_models, fR_Model(j) )
    #np.append(fR_models, fR_Model(j,fR_rad[j][:]))





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