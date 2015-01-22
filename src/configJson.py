'''
Created on Dec 10, 2014

@author: huongnvo
'''

import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import json

class configJson(Tkinter.Tk):  
    def __init__(self, parent, cc):
        Tkinter.Tk.__init__(self, parent)
        self.grid()
        self.parent = parent
        self.ee = cc
        self.dd = {}
        self.errorBox = None
        self.configJsonFile()
            
    def configJsonFile(self):     
        self.depthEntry      = Entry(self)
        self.adChoiceEntry   = Entry(self)
        self.siteChoiceEntry = Entry(self)
        self.z_resEntry      = Entry(self)
        self.stepsEntry      = Entry(self)
        self.userDataEntry   = Entry(self)
        self.gravityEntry    = Entry(self)
        self.thermalEntry    = Entry(self)
        self.pressureEntry   = Entry(self)
        self.conzoneEntry    = Entry(self)
        self.rho0Entry       = Entry(self)
        self.tGradEntry      = Entry(self)
        self.diffuEntry      = Entry(self)
        self.gasChoiceEntry  = Entry(self)
        self.runTypeEntry    = Entry(self)
          
        self.errorBox = Text(self, height = 30, width = 40)
        self.generateJson = Tkinter.Button(self, text = "Generate config file", height = 1, width = 20, command = self.generate)
        
        Label(self, text="Depth:").grid(row = 0, sticky = W)
        Label(self, text="Advection method: ").grid(row = 1, sticky = W)
        Label(self, text="Sitechoice: ").grid(row = 2, sticky = W)
        Label(self, text="Z-resolution: ").grid(row = 3, sticky = W)
        Label(self, text="Stepsize: ").grid(row = 4, sticky = W)
        Label(self, text="User's data: ").grid(row = 5, sticky = W)
        Label(self, text="Gravity: ").grid(row = 6, sticky = W)
        Label(self, text="Thermal: ").grid(row = 7, sticky = W)
        Label(self, text="Pressure: ").grid(row = 8, sticky = W)
        Label(self, text="ConZone depth: ").grid(row = 9, sticky = W)
        Label(self, text="Rho0: ").grid(row = 10, sticky = W)
        Label(self, text="T-grad: ").grid(row = 11, sticky = W)
        Label(self, text="Diffusion: ").grid(row = 12, sticky = W)
        Label(self, text="Gas choice: ").grid(row = 13, sticky = W)
        Label(self, text="Run type: ").grid(row = 14, sticky = W)
    
        self.depthEntry.grid(row = 0, column = 1)
        self.adChoiceEntry.grid(row = 1, column = 1)
        self.siteChoiceEntry.grid(row = 2, column =1)
        self.z_resEntry.grid(row = 3, column = 1)
        self.stepsEntry.grid(row = 4, column = 1)
        self.userDataEntry.grid(row = 5, column = 1)
        self.gravityEntry.grid(row = 6, column = 1)
        self.thermalEntry.grid(row = 7, column = 1)
        self.pressureEntry.grid(row = 8, column = 1)
        self.conzoneEntry.grid(row = 9, column = 1)
        self.rho0Entry.grid(row = 10, column = 1)
        self.tGradEntry.grid(row = 11, column = 1)
        self.diffuEntry.grid(row = 12, column = 1)
        self.gasChoiceEntry.grid(row = 13, column = 1)
        self.runTypeEntry.grid(row = 14, column = 1)
        self.errorBox.grid(row = 0, column = 2, rowspan = 15, pady = 5, padx = 5)
        self.generateJson.grid(row = 15, column = 0, columnspan = 13 )
        
        if not bool(self.ee):
            self.errorBox.insert(Tkinter.END, "Config file is currently empty\n")
        else:
            self.errorBox.insert(Tkinter.END, self.ee)      
        
    def generate(self):
        self.errorBox.insert(Tkinter.END, "Config file is currently being generated")
        
        depth       = float(self.depthEntry.get())
        variable1 = {"depth": depth}
        adchoice    = self.adChoiceEntry.get()
        variable2 = {"ad_method": adchoice}
        sitechoice  = self.siteChoiceEntry.get()
        variable3 = {"sitechoice": sitechoice}
        z_res       = float(self.z_resEntry.get())
        variable4 = {"z_res": z_res}
        steps       = float(self.stepsEntry.get())
        variable5 = {"steps": steps}
        userdata    = self.str2bool(self.userDataEntry.get())
        variable6 = {"userdata": userdata}
        gravity     = self.gravityEntry.get()
        variable7 = {"gravity": gravity}
        thermal     = self.thermalEntry.get()
        variable8 = {"thermal": thermal}
        pressure    = float(self.pressureEntry.get())
        variable9 = {"pressure": pressure}
        conzone     = float(self.conzoneEntry.get())
        variablea = {"conzone": conzone}
        rho0        = float(self.rho0Entry.get())
        variableb = {"rho0": rho0}
        tgrad       = float(self.tGradEntry.get())
        variablec = {"tgrad": tgrad}
        diffu       = self.diffuEntry.get()
        variabled = {"diffu": diffu}
        gaschoice   = self.gasChoiceEntry.get()
        variablee = {"gaschoice": gaschoice}
        runtype     = self.runTypeEntry.get()
        variablef = {"runtype": runtype}
        
        self.setConfig()(self.dd, variable1)
        self.setConfig()(self.dd, variable2)
        self.setConfig()(self.dd, variable3)
        self.setConfig()(self.dd, variable4)
        self.setConfig()(self.dd, variable5)
        self.setConfig()(self.dd, variable6)
        self.setConfig()(self.dd, variable7)
        self.setConfig()(self.dd, variable8)
        self.setConfig()(self.dd, variable9)
        self.setConfig()(self.dd, variablea)
        self.setConfig()(self.dd, variableb)
        self.setConfig()(self.dd, variablec)
        self.setConfig()(self.dd, variabled)
        self.setConfig()(self.dd, variablee)
        self.setConfig()(self.dd, variablef)
            
        self.save(self.dd)
        self.destroy()
    
    def setConfig(self):
        def field(dd, dic):
            if dic[dic.keys()[0]] != None:
                dd[dic.keys()[0]] = dic[dic.keys()[0]]
            elif self.ee[dic.keys()[0]] != None:
                dd[dic.keys()[0]] = self.ee[dic.keys()[0]]
            else:
                return
        return field
      
    def save(self, dd):
        self.f = tkFileDialog.asksaveasfile(mode = 'w', defaultextension = ".json")
        if self.f is None: # asksaveasfile return `None` if dialog closed with "cancel".
            return
        with open(str((self.f).name), 'w') as outfile:
            json.dump(dd, outfile, indent = 4, ensure_ascii = False)
            self.ran = True
        return outfile
    
    def str2bool(self, v):
        return v.lower() in ("yes", "true", "t", "1")
  
    def getConfig(self):
        return self.dd
    
    def getPath(self):
        return self.f
    
    def didRun(self):
        return self.ran
    
    def quit(self):
        global root
        root.quit()
    