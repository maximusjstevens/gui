'''
Created on Dec 10, 2014

@author: huongnvo
'''
import Tkinter
import tkFileDialog
from Tkinter import Menu
from Tkinter import *
import os
import json
import configJson
import firnAirModel

class airModel(Tkinter.Tk):

    def __init__(self, parent):
        Tkinter.Tk.__init__(self, parent)
        self.grid()
        self.parent = parent
        self.cc = {}
        self.initialize()   
        self.configRan = False
        self.dataRan = False
        self.saveRan = False
        
    def initialize(self):
        menubar = Menu(self)
        self.config(menu = menubar)
        filemenu = Menu(menubar, tearoff = 0)
        menubar.add_cascade(label = 'File', menu = filemenu)
        filemenu.add_command(label = "Exit", command = quit)
        filemenu.add_command(label = 'Manual', command = self.manual)
        
        jsonButton   = Tkinter.Button(text = "Select a config file", height = 1, width = 20, command = self.openConfig)    
        configButton = Tkinter.Button(text = "Make a config file", height = 1, width = 20, command = self.configJson)
        openButton   = Tkinter.Button(text = "Select a data directory", height = 1, width = 20, command = self.openData)
        saveButton   = Tkinter.Button(text = "Select a save directory", height = 1, width = 20, command = self.selectDir)
        self.runButton    = Tkinter.Button(text = "Run model", height = 1, width = 45, command = self.run, state = DISABLED)
        self.plotButton   = Tkinter.Button(text = "Generate graphs", height = 1, width = 45, command = self.display, state = DISABLED)
        self.loadButton   = Tkinter.Button(text = "Display results in larger window", height = 1, width = 45, command = self.display, state = DISABLED)
        clearButton  = Tkinter.Button(text = "Clear console", height = 1, width = 55, command = self.clear)
        
        self.text = Text(self.parent)
        self.text.insert(Tkinter.END, "Welcome to the Community Firn Model\n")
        self.text.insert(Tkinter.END, "\n")
        self.text.insert(Tkinter.END, "To get started, choose a data directory and a results directory\n")
        self.text.insert(Tkinter.END, "also choose a .json file or configure your own\n")
        self.text.insert(Tkinter.END, "\n")
    
        self.result = Text(self.parent, height = 15, width = 48)
          
        jsonButton.grid(row = 0, column = 0, padx = 5)
        configButton.grid(row = 0, column = 1, padx = 5) 
        openButton.grid(row = 1, column = 0, padx = 5)  
        saveButton.grid(row = 1, column = 1, padx = 5)
        
        self.runButton.grid(row = 2, column = 0, columnspan = 2)
        self.plotButton.grid(row = 3, column = 0, columnspan = 2)
        self.loadButton.grid(row = 4, column = 0, columnspan = 2)
        self.text.grid(row = 0, column = 2, rowspan = 35, pady = 5, padx = 10)
        self.result.grid(row = 5, column = 0, columnspan = 2, rowspan = 5, pady = 5, padx = 5)
        clearButton.grid(row = 36, column = 2, pady = 5)
    
    def configJson(self):
        config = configJson.configJson(None, self.cc)
        config.title("Configure your .json file")
        if (config.didRun()):
            self.cc = config.getConfig()
            self.f = config.getPath()
            self.write("Your .json file has been configured. The current values are:\n")
            self.write(self.cc)
        self.write("Finishing up\n")
        self.configRan = True
        self.checkState()
        
    def openConfig(self):
        file_opt = options = {} 
        options['filetypes'] = [('all files', '.*'), ('text files', '.txt'), ('json files', '.json')]
        options['title'] = 'Choose a config file'
        config = tkFileDialog.askopenfilename(**file_opt)
        if config is "":
            return
        with open(config) as f:
            json_data = f.read()
            self.cc = json.loads(json_data)
            self.write(config + " has been opened\n")
            self.write("Fields in .json file\n")
            self.write(self.cc)
            self.write("\n")
            self.write("\n")
        self.configRan = True
        self.checkState()
        
    def openData(self):
        self.dataDir = tkFileDialog.askdirectory()
        if self.dataDir is "":
            return
        self.write(self.dataDir)
        files = os.listdir(self.dataDir)
        self.write(self.dataDir + " has been opened\n")
        self.write("The following files have been loaded:\n")
        for dataFile in files:
            self.write("    " + dataFile + "\n")
        self.write("\n")
        self.dataRan = True
        self.checkState()
        
    def selectDir(self):
        self.resultDir = tkFileDialog.askdirectory()
        if self.resultDir is "":
            self.write("You did not choose a directory")
        else:
            self.write("Results will be saved to: " + self.resultDir)
        self.saveRan = True
        self.checkState()
    
    def getDir(self):
        return self.resultDir
        
    def write(self, txt):
        self.text.insert(Tkinter.END, str(txt))

    def clear(self):
        self.text.delete(1.0, Tkinter.END)
    
    def display(self):
        print "display"
        
    def run(self):
        print(str(self.dataDir))
        print(str(self.resultDir))
        firnAir = firnAirModel.firnAirModel(self.cc, self.dataDir, self.resultDir)
    
    def saveData(self):
        print "saving data"
    
    def manual(self):
        print "print manual"
        
    def checkState(self):
        if self.dataRan and self.saveRan and self.configRan:
            self.runButton['state'] = 'normal'
            self.loadButton['state'] = 'normal'
            self.plotButton['state'] = 'normal'
     
if __name__ == "__main__":
    air = airModel()
