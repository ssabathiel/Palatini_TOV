import sys
from PyQt5.QtWidgets import QDialog, QApplication, QPushButton, QVBoxLayout, QCheckBox, QHBoxLayout
from PyQt5.QtCore import Qt

#import PyQt4
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtWidgets import QTreeWidget, QTreeWidgetItem, QApplication, QWidget
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.image as image
matplotlib.rcParams['backend'] = "TkAgg"
import random
import glob
from Models import GR_Model
from Models import fR_Model
from Models import fRQ_Model
import subprocess
from itertools import cycle


lines = ["--","-.",":"]
linecycler = cycle(lines)


sil_fR=1
sil_fRQ=1
sil_GR=1
diff=0
profile = 0

plot_title = 'f(R),f(R,Q),GR comparison for Mass and Radius'


###############
## Get data from files to arrays
###############


fR_Rp = glob.glob('Results/TOV_output_fR_*')[0][25:32] #35]
fR_Rp +=   '^{' + glob.glob('Results/TOV_output_fR_*')[0][33:35] + '}'
print("Rp=", fR_Rp)

fRQ_Rp = glob.glob('Results/TOV_output_fRQ_*')[0][26:33] #35]
fRQ_Rp +=   '^{' + glob.glob('Results/TOV_output_fRQ_*')[0][34:36] + '}'
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


fR_Modellist = [fR_Model (count) for count in range(len(glob.glob('Results/TOV_output_fR_*')))]
fRQ_Modellist = [fRQ_Model (count) for count in range(len(glob.glob('Results/TOV_output_fRQ_*')))]
GR_Modellist = [GR_Model (count) for count in range(len(glob.glob('Results/TOV_output_GR_*')))]


for i in range(len(glob.glob('Results/TOV_output_fR_*'))):
    filename_length = len(glob.glob('Results/TOV_output_fR_*')[i])
    fR_Modellist[i].EOS_type = glob.glob('Results/TOV_output_fR_*')[i][-3:filename_length]

    fR_Modellist[i].ID =  glob.glob('Results/TOV_output_fR_*')[i][25:32]
    fR_Modellist[i].ID += '^{' + glob.glob('Results/TOV_output_fR_*')[i][33:35] + '}'
    with open(glob.glob('Results/TOV_output_fR_*')[i]) as f:
        lines = f.readlines()
        fR_Modellist[i].rad =  [line.split()[2] for line in lines]
        fR_Modellist[i].mass = [line.split()[1] for line in lines]
        fR_Modellist[i].rho = [line.split()[0] for line in lines]

for i in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
    filename_length = len(glob.glob('Results/TOV_output_fRQ_*')[i])
    fRQ_Modellist[i].EOS_type = glob.glob('Results/TOV_output_fRQ_*')[i][-3:filename_length]
    fRQ_Modellist[i].ID = glob.glob('Results/TOV_output_fRQ_*')[i][26:33]  # 35]
    fRQ_Modellist[i].ID += '^{' + glob.glob('Results/TOV_output_fRQ_*')[i][34:36] + '}'
    fRQ_Modellist[i].ID += ' R_q= '
    fRQ_Modellist[i].ID += glob.glob('Results/TOV_output_fRQ_*')[i][41:48]
    fRQ_Modellist[i].ID += '{' + glob.glob('Results/TOV_output_fRQ_*')[i][48:50] + '}'
    with open(glob.glob('Results/TOV_output_fRQ_*')[i]) as f:
        lines = f.readlines()
        fRQ_Modellist[i].rad = [line.split()[2] for line in lines]
        fRQ_Modellist[i].mass = [line.split()[1] for line in lines]
        fRQ_Modellist[i].rho = [line.split()[0] for line in lines]

for i in range(len(glob.glob('Results/TOV_output_GR_*'))):
    filename_length = len(glob.glob('Results/TOV_output_GR_*')[i])
    GR_Modellist[i].ID =  glob.glob('Results/TOV_output_GR_*')[i][-3:filename_length ]
    GR_Modellist[i].EOS_type =  glob.glob('Results/TOV_output_GR_*')[i][-3:filename_length ]
    #GR_Modellist[i].ID += '{' + glob.glob('Results/TOV_output_GR_*')[i][32:35] + '}'
    with open(glob.glob('Results/TOV_output_GR_*')[i]) as f:
        lines = f.readlines()
        GR_Modellist[i].rad = [line.split()[2] for line in lines]
        GR_Modellist[i].mass = [line.split()[1] for line in lines]
        GR_Modellist[i].rho = [line.split()[0] for line in lines]


for j in range(len(glob.glob('Results/TOV_output_fR_*'))):
    for i in range(len(fR_Modellist[j].rad)):
        fR_Modellist[j].rad[i] = float(fR_Modellist[j].rad[i])/100000.
        fR_Modellist[j].rho[i] = np.log10(float(fR_Modellist[j].rho[i]))

for j in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
    for i in range(len(fRQ_Modellist[j].rad)):
        fRQ_Modellist[j].rad[i] = float(fRQ_Modellist[j].rad[i])/100000.
        fRQ_Modellist[j].rho[i] = np.log10(float(fRQ_Modellist[j].rho[i]))

for j in range(len(glob.glob('Results/TOV_output_GR_*'))):
    for i in range(len(GR_Modellist[j].rad)):
        GR_Modellist[j].rad[i] = float(GR_Modellist[j].rad[i])/100000.
        GR_Modellist[j].rho[i] = np.log10(float(GR_Modellist[j].rho[i]))




class Window(QDialog):
    def __init__(self, parent=None):
        super(Window, self).__init__(parent)



        self.plot_title = 'f(R),f(R,Q),GR comparison for Mass and Radius'
        self.xlabel = "$R[km]$"
        self.ylabel = r'$M[M_\odot]$'
        self.xlim = (8.0, 20.0)

        tree = QTreeWidget()
        headerItem = QTreeWidgetItem()
        item = QTreeWidgetItem()

        tree_fRQ = QTreeWidget()
        headerItem_fRQ = QTreeWidgetItem()
        item_fRQ = QTreeWidgetItem()

        tree_GR = QTreeWidget()
        headerItem_GR = QTreeWidgetItem()
        item_GR = QTreeWidgetItem()

        tree_quantities = QTreeWidget()
        headerItem_quantities = QTreeWidgetItem()
        item_quantities = QTreeWidgetItem()

        self.quant_branch = QTreeWidgetItem(tree_quantities)
        self.quant_branch.setText(0, "Quantities")
        self.quant_branch.setFlags(self.quant_branch.flags() )
        # self.GR_branch.setCheckState(0, Qt.Unchecked)

        child = QTreeWidgetItem(self.quant_branch)
        child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
        child.setText(0, 'R-M')
        child.setCheckState(0, Qt.Checked)
        child.emitDataChanged()  # .connect(self.plot)

        child = QTreeWidgetItem(self.quant_branch)
        child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
        child.setText(0, 'rho-M')
        child.setCheckState(0, Qt.Unchecked)
        child.emitDataChanged()

        self.GR_branch = QTreeWidgetItem(tree_GR)
        self.GR_branch.setText(0, "GR")
        self.GR_branch.setFlags(self.GR_branch.flags() | Qt.ItemIsTristate  | Qt.ItemIsUserCheckable)
        #self.GR_branch.setCheckState(0, Qt.Unchecked)

        for x in range(len(glob.glob('Results/TOV_output_GR_*'))):
            child = QTreeWidgetItem(self.GR_branch)
            child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
            checktext =  GR_Modellist[x].EOS_type
            child.setText(0,checktext)
            child.setCheckState(0, Qt.Unchecked)
            child.emitDataChanged()# .connect(self.plot)

        self.fR_branch = QTreeWidgetItem(tree)
        self.fR_branch.setText(0,"fR")
        self.fR_branch.setFlags(self.fR_branch.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)

        for x in range(len(glob.glob('Results/TOV_output_fR_*'))):
            child = QTreeWidgetItem(self.fR_branch)
            child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
            checktext = fR_Modellist[x].ID + fR_Modellist[x].EOS_type
            child.setText(0,checktext)
            child.setCheckState(0, Qt.Unchecked)
            child.emitDataChanged()# .connect(self.plot)

        #if(fRQ_branch.child(1).isSelected==1){print("hey hey")}

        self.fRQ_branch = QTreeWidgetItem(tree_fRQ)
        self.fRQ_branch.setText(0, "fRQ")
        self.fRQ_branch.setFlags(self.fRQ_branch.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)

        for x in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
            child = QTreeWidgetItem(self.fRQ_branch)
            child.setFlags(child.flags() | Qt.ItemIsUserCheckable)
            checktext = fRQ_Modellist[x].ID + fRQ_Modellist[x].EOS_type
            child.setText(0, checktext)
            child.setCheckState(0, Qt.Unchecked)
        print("checked?", self.fR_branch.child(0).isSelected())
        tree.itemChanged.connect(self.plot)
        tree_fRQ.itemChanged.connect(self.plot)
        tree_GR.itemChanged.connect(self.plot)
        tree_quantities.itemChanged.connect(self.quant_change)
        tree_quantities.itemChanged.connect(self.plot)



        #tree.show()

        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # Just some button connected to `plot` method
        #self.button = QPushButton('Plot')
        #self.button.clicked.connect(self.plot)

        # set the layout
        self.layout = QVBoxLayout()

        self.sublayout1 = QVBoxLayout()
        self.sublayout2 = QHBoxLayout()
        plotBox = QHBoxLayout()
        plotBox.addLayout(self.sublayout1, 1)
        plotBox.addLayout(self.sublayout2, 4)
        bottomlayout = QHBoxLayout()
        #self.layout.addLayout()


        self.layout.addWidget(self.toolbar)
        self.layout.addWidget(self.canvas,1)
        #self.layout.addWidget(self.button)
        self.setLayout(self.layout)

        #layout = QHBoxLayout()
        self.btn_rhoM = QCheckBox("rho-M")
        self.btn_RM = QCheckBox("R-M")
        self.btn_GR = QCheckBox("GR")
        self.btn_fR =  QCheckBox("fR")
        self.btn_fRQ = QCheckBox("fRQ")

        self.btn_fRs = []

        self.btn_rhoM.setChecked(False)
        self.btn_RM.setChecked(True)
        self.btn_GR.setChecked(False)
        self.btn_fR.setChecked(False)
        self.btn_fRQ.setChecked(False)

        #self.b1.stateChanged.connect(lambda: self.btnstate(self.b1))
        #self.btn_rhoM.stateChanged.connect(self.quant_change)
        #self.btn_RM.stateChanged.connect(self.quant_change)

        self.btn_rhoM.stateChanged.connect(self.plot)
        self.btn_RM.stateChanged.connect(self.plot)
        self.btn_GR.stateChanged.connect(self.plot)
        self.btn_fR.stateChanged.connect(self.plot)
        self.btn_fRQ.stateChanged.connect(self.plot)

        #self.layout.addWidget(self.btn_GR)
        #self.layout.addWidget(self.btn_fR)
        #self.layout.addWidget(self.btn_fRQ)
        self.sublayout2.addWidget(tree_quantities)
        self.sublayout2.addWidget(tree_GR)
        self.sublayout2.addWidget(tree)
        self.sublayout2.addWidget(tree_fRQ)
        self.layout.addLayout(plotBox)

    def quant_change(self):
        if(self.quant_branch.child(0).checkState(0)):
            self.plot_title = 'f(R),f(R,Q),GR comparison for Radius and Mass of NS'
            self.xlabel = "R[km]"
            self.ylabel = r'M[$M_\odot$]'
            self.xlim = (8.0, 20.0)


        if (self.quant_branch.child(1).checkState(0)):
            self.plot_title = 'f(R),f(R,Q),GR comparison for Central Density and Mass of NS'
            self.xlabel = r'$\log\rho_c[\frac{g}{cm^3}]$'
            self.ylabel = r'M[$M_\odot$]'
            self.xlim = (13.5, 17.0)


    def plot(self):
        fR_label = 'f(R)\n  $(R_p=$' + r'$%s$' % fR_Rp + r'$)$'
        fRQ_label = 'f(R,Q) \n  $(R_p=$' + r'$%s$' % fRQ_Rp + '$)$\n  ' + r'$(R_q=$' + r'$%s$' % fRQ_Rq + '$)$'

        if self.btn_fR.isChecked() == True:
            #self.btn_fRs.append(QCheckBox("yep"))
            print("self.fR_branch.child(i).isSelected()", self.fR_branch.child(i).isSelected())
            self.layout.addWidget(self.btn_fRs[0])
            plt.plot(x, y,  label=r'${}$'.format(fR_label), color='blue')
            plt.title('fR/GR comparison for Mass and Radius for $R_p=2.1*10^{80}$')
            plt.legend(fontsize=20)

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

        for i in range(len(glob.glob('Results/TOV_output_GR_*'))):
            if self.GR_branch.child(i).checkState(0):
                labelly = 'GR ' + GR_Modellist[i].EOS_type
                if (self.quant_branch.child(0).checkState(0)):
                    ax.plot(GR_Modellist[i].rad, GR_Modellist[i].mass, label=labelly) #, next(linecycler)
                if (self.quant_branch.child(1).checkState(0)):
                    ax.plot(GR_Modellist[i].rho, GR_Modellist[i].mass, label=labelly) #, next(linecycler)



        for i in range(len(glob.glob('Results/TOV_output_fR_*'))):
            if self.fR_branch.child(i).checkState(0):
                labelly = 'f(R) ' + fR_Modellist[i].EOS_type + ": R_p= " + fR_Modellist[i].ID
                if (self.quant_branch.child(0).checkState(0)):
                    ax.plot(fR_Modellist[i].rad,fR_Modellist[i].mass,next(linecycler),label=r'${}$'.format(labelly))
                if (self.quant_branch.child(1).checkState(0)):
                    ax.plot(fR_Modellist[i].rho, fR_Modellist[i].mass,next(linecycler), label=r'${}$'.format(labelly))
                #print(self.fR_branch.child(i).isSelected())
                #ax.plot(x5, y5, label='GR', color='black')


        for i in range(len(glob.glob('Results/TOV_output_fRQ_*'))):
            if self.fRQ_branch.child(i).checkState(0):
                labelly = 'f(R,Q) ' + fRQ_Modellist[i].EOS_type + ": R_p= " + fRQ_Modellist[i].ID
                if (self.quant_branch.child(0).checkState(0)):
                    ax.plot(fRQ_Modellist[i].rad,fRQ_Modellist[i].mass,next(linecycler),label=r'${}$'.format(labelly))
                if (self.quant_branch.child(1).checkState(0)):
                    ax.plot(fRQ_Modellist[i].rho, fRQ_Modellist[i].mass, next(linecycler), label=r'${}$'.format(labelly))

                    #ax.plot(x5, y5, label='GR', color='black')
                #print(fRQ_Modellist[i].rad, " ", fRQ_Modellist[i].mass)


        ax.legend(fontsize=15)
        if self.btn_fRQ.isChecked() == True:
            ax.plot(x2, y2, label=fRQ_label, color='blue')
            ax.legend(fontsize=20)



        #ax.set_title(self.plot_title)
        ax.set_xlim(self.xlim)
        ax.set_xlabel(self.xlabel,fontsize=20)
        ax.set_ylabel(self.ylabel,fontsize=20)
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