import numpy as np
import matplotlib



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
        self.ID = ""
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
