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
        self.resize(900,600)

        # plotWidget1: Dit versus phis distribution 
        self.plotWidget1 = pg.PlotWidget()
        self.plotWidget1.setBackground("w")
        self.plotWidget1.setTitle("Dit distribution", color="b", size="10pt")
        self.plotWidget1.addLegend(offset=(-10, 10))
        styles={'color':'b'}
        self.plotWidget1.setLabel("left", "log10 Dit", "/cm2/eV", **styles)
        self.plotWidget1.setLabel("bottom", "E-Ei", "eV", **styles)
        self.plotWidget1.showGrid(x=True, y=True)
        self.plotWidget1.setRange(xRange=(-2,2.0), yRange=(10, 13), padding=0)
        
        layout1=QtWidgets.QVBoxLayout()
        layout1.setContentsMargins(0, 0, 0, 0)
        layout1.addWidget(self.plotWidget1)
        self.graphicsView_Dit.setLayout(layout1)
        
        # plotWidget2: Charge versus phis 
        self.plotWidget2 = pg.PlotWidget()
        self.plotWidget2.setBackground("w")
        self.plotWidget2.setTitle("Charge", color="b", size="10pt")
        self.plotWidget2.addLegend(offset=(-10, 10))
        styles={'color':'b'}
        self.plotWidget2.setLabel("left", "log10(|Q|)", "C/cm2", **styles)
        self.plotWidget2.setLabel("bottom", "phis", "eV", **styles)
        self.plotWidget2.showGrid(x=True, y=True)
        self.plotWidget2.setRange(xRange=(-2.0, 2.0), yRange=(-9, -4), padding=0)

        layout2=QtWidgets.QVBoxLayout()
        layout2.setContentsMargins(0, 0, 0, 0)
        layout2.addWidget(self.plotWidget2)
        self.graphicsView_Q.setLayout(layout2)

        # plotWidget3: C-V
        self.plotWidget3 = pg.PlotWidget()
        self.plotWidget3.setBackground("w")
        self.plotWidget3.setTitle("capacitance", color="b", size="10pt")
        self.plotWidget3.addLegend(offset=(-80,50))
        styles={'color':'b'}
        self.plotWidget3.setLabel("left", "C", "F/cm2", **styles)
        self.plotWidget3.setLabel("bottom", "VG", "V", **styles)
        self.plotWidget3.showGrid(x=True, y=True)
        self.plotWidget3.setRange(xRange=(-10.0,10.0), yRange=(0, 50e-9), padding=0)
        
        layout3=QtWidgets.QVBoxLayout()
        layout3.setContentsMargins(0, 0, 0, 0)
        layout3.addWidget(self.plotWidget3)
        self.graphicsView_Capacitance.setLayout(layout3)

        self.pushButton_calculate.clicked.connect(self.calculate)
        self.pushButton_clear.clicked.connect(self.clearGraph)

        self.calculate()

    def calculate(self):
        beta=QE/KB/TEMP
        self.NA=(self.doubleSpinBox_NA_Mantissa.value())*10**float((self.spinBox_NA_Exponent.value()))
        self.NI=(self.doubleSpinBox_NI_Mantissa.value())*10**float((self.spinBox_NI_Exponent.value()))
        self.PsiB=math.log(math.fabs(self.NA)/self.NI)/beta # PsiB: Fermi level of substrate from Ei
        if self.NA<0:
            self.PsiB=-self.PsiB # n-type PsiB < 0 (Sze style)
        self.label_PsiB.setText(str(self.PsiB))
        self.phis1=self.doubleSpinBox_phis1.value()
        self.phis2=self.doubleSpinBox_phis2.value()
        self.phisStep=(self.doubleSpinBox_phisStep.value())*1e-3 # (meV)
        self.tox=(self.doubleSpinBox_DielectricThickness.value())*1e-7 # (nm)
        self.epsd=self.doubleSpinBox_DielectricPermittivity.value()*EPS0

        self.Dit_Ecnl=self.doubleSpinBox_Dit_Ecnl.value()
        self.Dit_Dit0=(self.doubleSpinBox_Dit_Dit0_Mantissa.value())*10**float((self.spinBox_Dit_Dit0_Exponent.value()))
        self.Dit_a=self.doubleSpinBox_Dit_a.value()
        self.Dit_c=self.doubleSpinBox_Dit_c.value()
        
        self.draw_Dit()
        self.draw_Q_CV()
    
    def draw_Q_CV(self):
        beta=QE/KB/TEMP
        phis1=self.phis1
        phis2=self.phis2
        phisStep=self.phisStep
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
#        Qs=[]
        Qslog_posi=[]
        Qslog_nega=[]
        CD=[]
        Ctotal=[]
        VG=[]
        VGDit=[]
#        phis = np.arange(phis1, phis2, phisStep)
        phis = frange(phis1, phis2, phisStep)
        Cox=self.epsd/self.tox
#        print(f'phis1 {phis1:.3g}, phis2 {phis2:.3g} phisStep {phisStep:.3g} Cox {Cox:.3g}')
        for phis0 in phis:
            Qs0=-(math.sqrt(2.0)*EPSS/beta/LDpp0 * Ffunction(beta*phis0, np0/pp0))
            if phis0 < 0:
                Qs0=-Qs0
            if Qs0>0:
                Qslog_nega.append(math.nan)
                try:
                    Qslog_posi.append(math.log10(Qs0))
                except:
                    Qslog_posi.append(math.nan)
            else:
                Qslog_posi.append(math.nan)
                try:
                    Qslog_nega.append(math.log10(-Qs0))
                except:
                    Qslog_nega.append(math.nan)

            try:
                CD0=( EPSS/math.sqrt(2.0)/LDpp0 *math.fabs(1.0-math.exp(-beta*phis0)+(np0/pp0)*(math.exp(beta*phis0)-1.0)) / Ffunction(beta*phis0, np0/pp0) )
            except ZeroDivisionError:
                if self.NA>0:
                    CD0=( EPSS/LDpp0 )
                else:
                    CD0=( EPSS/LDnp0 )
            CD.append(CD0)
            Ctotal.append(1.0/ (1.0/Cox + 1.0/CD0) )

            Qit0=self.interfaceCharge(phis0)
            VG.append((phis0  - (Qs0)/Cox))
            VGDit.append((phis0 - (Qs0+Qit0)/Cox))

        # Plot Qs - phis
        pen=pg.mkPen('b', width=2, style=QtCore.Qt.SolidLine)
        self.plotWidget2.plot(x=phis, y=Qslog_posi, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="posi |Qs|", ignoreBounds=True)
        pen=pg.mkPen('b', width=2, style=QtCore.Qt.DashLine)
        self.plotWidget2.plot(x=phis, y=Qslog_nega, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="nega |Qs|")

        # Plot C-V
        VG_clip=np.clip(VG, -10, 10) # clip is recommended for stable drawing
        pen=pg.mkPen(color="b")
        self.plotWidget3.plot(x=VG_clip, y=Ctotal, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="Ctotal w/o Dit", ignoreBounds=True)

        # Plot C-V with Dit
        VGDit_clip=np.clip(VGDit, -10, 10) # clipping is recommended for stable drawing
        pen=pg.mkPen(color="r")
        self.plotWidget3.plot(x=VGDit_clip, y=Ctotal, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="Ctotal w/ Dit", ignoreBounds=True)

        Qitlog_posi=[]
        Qitlog_nega=[]
        for phis0 in phis:
            try:
                Qitlog_posi.append(math.log10(self.interfaceCharge(phis0)))
            except:
                Qitlog_posi.append(math.nan)
            try:
                Qitlog_nega.append(math.log10(-self.interfaceCharge(phis0)))
            except:
                Qitlog_nega.append(math.nan)

        # Plot Qit-V
        pen=pg.mkPen('r', width=2, style=QtCore.Qt.SolidLine)
        self.plotWidget2.plot(x=phis, y=Qitlog_posi, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="posi |Qit|")
        pen=pg.mkPen('m', width=2, style=QtCore.Qt.DashLine)
        self.plotWidget2.plot(x=phis, y=Qitlog_nega, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="nega |Qit|")
        
        
    def clearGraph(self):
        self.plotWidget1.clear()
        self.plotWidget2.clear()
        self.plotWidget3.clear()
        pg.QtGui.QGuiApplication.processEvents()

    def draw_Dit(self):
        phis = np.arange(-1.0, 1.0, 0.02)
        E_Ei=[]
        Dit=[]
        Ditlog_posi=[]
        Ditlog_nega=[]
        for phis0 in phis:
            E_Ei.append(phis0)
            Dit0=self.interfaceStateDensity(phis0)
            Dit.append(Dit0)
            if Dit0>0:
                Ditlog_nega.append(math.nan)
                try:
                    Ditlog_posi.append(math.log10(Dit0))
                except:
                    Ditlog_posi.append(math.nan)
            else:
                Ditlog_posi.append(math.nan)
                try:                    
                    Ditlog_nega.append(math.log10(-Dit0))
                except:
                    Ditlog_nega.append(math.nan)

        # Plot Dit
        pen=pg.mkPen('r', width=2, style=QtCore.Qt.SolidLine)
        self.plotWidget1.plot(x=E_Ei, y=Ditlog_posi, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="posi Dit")
        pen=pg.mkPen('m', width=2, style=QtCore.Qt.DashLine)
        self.plotWidget1.plot(x=E_Ei, y=Ditlog_nega, pen=pen, brush=pg.mkBrush("b"), size =7.5, antialias = True, name="nega Dit")

    def interfaceStateDensity(self, efs): # Dit (/cm2/eV)
        # efs: surface Fermi level from Ei
        # Dit model function: dit=Dit0*exp(c*|((e-ecnl)|**a)
        # Dit0: Dit minimun (/cm^2/eV)
        # a: factor (cm^2/eV)
        # c: shape factor
        # ecnl: E_(charge neutral level) from Ei (Ec side postivie)
        ecnl=self.Dit_Ecnl
        dit0=self.Dit_Dit0
        a=self.Dit_a
        c=self.Dit_c
        if (efs > ecnl):
            sign=-1.0 # accepter-like interface state
        else:
            sign=+1.0 # donnor-like interface state
            
        return sign*(dit0*math.exp(c*pow((math.fabs(efs-ecnl)),a))) # (/cm2/eV)
            
    def interfaceCharge(self, psis): # Qit C/cm2 */
        # psis surface potential
        # efs temporary surface Fermi level from Ei
        # ecnl from Ei (Ec side positive)
        # efsstep=1e-5; /* 0.01 mV step */
        # Qit;
        #  const double eNL0=Dit_eNL0;
        efs_destination=+psis-self.PsiB # from Ei
        ecnl=self.Dit_Ecnl
        dit0=self.Dit_Dit0
        a=self.Dit_a
        c=self.Dit_c

        Qit=0;
        efsstep=1e-3 # 0.01 mV step
        efs0=ecnl # integration from E_CNL (from Ei) to efs
        if (efs_destination > ecnl): # integration up or down?
            while (efs0 < (efs_destination)): # integration upward
                Qit=Qit+QE*self.interfaceStateDensity(efs0)*efsstep
                efs0+=efsstep
        else:
            while (efs0 > (efs_destination)): # integration downward
                Qit=Qit+QE*self.interfaceStateDensity(efs0)*efsstep
                efs0-=efsstep
        return Qit # (C/cm2)
        
def Ffunction(betaphi, npratio):
    try:
        return math.sqrt((math.exp(-betaphi) + betaphi -1.0) + npratio*((math.exp(betaphi) - betaphi -1.0)))
    except:
        return 0

def frange(start, stop, step):
    n = int(round((stop - start) / step))
    return [start + i*step for i in range(n+1)]

app=QtWidgets.QApplication(sys.argv)
main=MainWindow()
main.show()

app.exec_()
