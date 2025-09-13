import os, sys, math

from PyQt5 import QtGui, QtWidgets, uic
from PyQt5.QtWidgets import (
    QFileDialog, QDialog, QApplication,
    QPushButton, QMainWindow, QMessageBox
    )
    
from pyqtgraph.Qt import QtCore
#from pyqtgraph import LogAxisItem
import pyqtgraph as pg
import numpy as np


# Physical Constants
QE = 1.60218e-19 # (C) elementary charge
KB = 1.38066e-23 # (J/K) Plack's constant
TEMP = 300.0       # (K) temperature
EPS0 = 8.85418e-14 # (F/cm) permittivity in vacuum


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        uic.loadUi("CVCalculator_gui.ui", self)
        self.initUI()

    def initUI(self):
        self.resize(800,600)

        # setting for Dit versus phis distribution 
        self.graphicsView_Dit.setEnabled(True)
        self.graphicsView_Dit.setBackground("white")
        self.graphicsView_Dit.setTitle("Dit distribution", color="b", size="10pt")
        styles={'color':'b'}
        self.graphicsView_Dit.setLabel("left", "Dit", "/cm2/eV", **styles)
        self.graphicsView_Dit.setLabel("bottom", "phis", "eV", **styles)

        # setting for Charge versus phis 

        self.graphicsView_Q.setEnabled(True)
        self.graphicsView_Q.setBackground("white")
        self.graphicsView_Q.setTitle("Charge", color="b", size="10pt")
        self.graphicsView_Q.setLabel("left", "|Q|", "C/cm2", **styles)
        self.graphicsView_Q.setLabel("bottom", "phis", "eV", **styles)
        self.graphicsView_Q.setLogMode(True)
        self.graphicsView_Q.setLogMode(x=False, y=False)
#        log_y_axis=LogAxisItem(orientation='left')
#        self.graphicsView_Q.PlotWidget(axisItems={'left': log_y_axis})
        self.graphicsView_Q.showGrid(x=True, y=True)
        self.graphicsView_Q.setRange(xRange=(-2.0, 2.0), yRange=(-9, -4), padding=0)

        self.graphicsView_Capacitance.setEnabled(True)
        self.graphicsView_Capacitance.setBackground("w")
        self.graphicsView_Capacitance.setTitle("capacitance", color="b", size="10pt")
        self.graphicsView_Capacitance.setLabel("left", "C", "F", **styles)
        self.graphicsView_Capacitance.setLabel("bottom", "VG", "V", **styles)
        self.graphicsView_Capacitance.setRange(xRange=(-10.0,10.0), yRange=(0, 100e-9), padding=0)


        self.pushButton_calculate.clicked.connect(self.calculate)

        self.calculate()

    def calculate(self):
        self.NA=(self.doubleSpinBox_NA_Mantissa.value())*10**float((self.spinBox_NA_Exponent.value()))
        self.NI=(self.doubleSpinBox_NI_Mantissa.value())*10**float((self.spinBox_NI_Exponent.value()))
        self.phis1=self.doubleSpinBox_phis1.value()
        self.phis2=self.doubleSpinBox_phis2.value()
        self.tox=self.doubleSpinBox_DielectricThickness.value()*1e-7
        self.epsd=self.doubleSpinBox_DielectricPermittivity.value()*EPS0

#        self.graphicsView_Dit = pg.PlotWidget()
#        self.graphicsView_Dit.plot(hour,temp)
        self.draw_Dit()
        self.draw_Q_CV()
#        self.draw_CV()
    
    def draw_Dit(self):
        pen=pg.mkPen(color="b")
        hour0=[1,2,3,4,5]
        temp0=[33,35,28,36,30]
#        self.graphicsView_Dit.plot(x=hour, y=temp,pen,"+",symbolSize=30, symbolBrush="b")
        self.graphicsView_Dit.addItem(pg.PlotCurveItem(x=hour0, y=temp0, symbol="x", pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True))


    def draw_Q_CV(self):
        beta=QE/KB/TEMP
        phis1=self.phis1
        phis2=self.phis2
        EPSS=11.9*EPS0
        if self.NA>0:
            pp0=self.NA
            np0=self.NI*self.NI/pp0
        else:
            np0=math.fabs(self.NA)
            pp0=self.NI*self.NI/np0
        npratio=np0/pp0
        LDpp0=math.sqrt(EPSS/QE/pp0/beta)
        LDnp0=math.sqrt(EPSS/QE/np0/beta)
#        print(f"pp0 {pp0:.3g} np0 {np0:.3g} npratio {npratio:.3g}")
#        print(self.NA, LD, npratio)
        Qs=[]
        Qslog=[]
        CD=[]
        Ctotal=[]
        VG=[]
        phis = np.arange(phis1, phis2, 0.02)
        Cox=self.epsd/self.tox
        for phis0 in phis:
            Qs0=-(math.sqrt(2.0)*EPSS/beta/LDpp0 * Ffunction(beta*phis0, np0/pp0))
            if phis0 > 0:
                Qs0=-Qs0
            Qs.append(Qs0)
            if math.fabs(Qs0)>0:
                Qslog.append(math.log10(math.fabs(Qs0)))
            else:
                Qslog.append(math.log10(1e-30))

            if math.fabs(phis0)>1e-6:
                CD0=( EPSS/math.sqrt(2.0)/LDpp0 *math.fabs(1.0-math.exp(-beta*phis0)+(np0/pp0)*(math.exp(beta*phis0)-1.0)) / Ffunction(beta*phis0, np0/pp0) )
                CD.append(CD0)
            else:
                if self.NA>0:
                    CD.append( EPSS/LDpp0)
                else:
                    CD.append( EPSS/LDnp0)
#        print(f"Cox {Cox:.3g}")
            Ctotal.append(1.0/ (1.0/Cox+1.0/CD0) )
#        print(CD)

            VG.append((phis0  + (Qs0)/Cox))
#        VG=[Qs0 for phis0, Qs0 in zip(phis, Qs)]
        print(Qs)

        # Plot Qs - phis
        pen=pg.mkPen(color="b")
        self.graphicsView_Q.addItem(pg.PlotCurveItem(x=phis, y=Qslog, symbol="x", pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True))

        # Plot C-V
        pen=pg.mkPen(color="b")
        self.graphicsView_Capacitance.addItem(pg.PlotCurveItem(x=VG, y=Ctotal, symbol="x", pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True))
        

def Ffunction(betaphi, npratio):
    try:
        return math.sqrt((math.exp(-betaphi) + betaphi -1.0) + npratio*((math.exp(betaphi) - betaphi -1.0)))
    except:
        return 0

app=QtWidgets.QApplication(sys.argv)
main=MainWindow()
main.show()

app.exec_()
