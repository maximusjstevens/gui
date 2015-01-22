'''
Created on Oct 6, 2014

@author: huongnvo
'''

import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import firnModel
import airModel

class mainModel(Tkinter.Tk):
    def __init__(self, parent):
        Tkinter.Tk.__init__(self, parent)
        self.grid()
        self.parent = parent
        self.initialize()
    
    def initialize(self):
        menubar = Menu(self)
        self.config(menu = menubar)
        filemenu = Menu(menubar, tearoff = 0)
        menubar.add_cascade(label = "File", menu = filemenu)
        filemenu.add_command(label = "Exit", command = quit)
    
        firnModel  = Tkinter.Button(text = "Run the firn model", height = 1, width = 50, command = self.firnModel)
        airModel   = Tkinter.Button(text = "Run the air model", height = 1, width = 50, command = self.airModel)    
        firnAndAir = Tkinter.Button(text = "Run both the firn & air model", height = 1, width = 50) 
             
        self.text = Text(self.parent, height = 15, width = 75)
        self.text.insert(Tkinter.END, "Welcome to the Community Firn Model\n")
        self.text.insert(Tkinter.END, "To get started, choose the model you would like to run")        
        
        firnModel.grid(row = 1, column = 0, padx = 5)
        airModel.grid(row = 2, column = 0, padx = 5)
        firnAndAir.grid(row = 3, column = 0, padx = 5)     
        self.text.grid(row = 0, column = 0, pady = 5, padx = 10) 
    
    def firnModel(self):
        self.destroy()
        firn = firnModel.firnModel(None)
        firn.title("Firn Model")
    
    def airModel(self):
        self.destroy()
        air = airModel.airModel(None)
        air.title("Air Model")
        
if __name__ == "__main__":
    app = mainModel(None)
    app.title("Community Firn Model")
    app.mainloop()
    
    
    

        
