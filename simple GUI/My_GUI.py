# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 14:57:42 2020

@author: hohle
"""

from tkinter import *
from tkinter.filedialog import askopenfilename
from MyModule import my_data_analysis as my_data_analysis
from tkinter import ttk

my_window = Tk()#creates a window
my_window.title('My first GUI')
my_window.geometry('850x150') #size in pixel


###############################################################################
#button for file selection
file_choose_button = Button(my_window,text = 'select a file',\
                            command = lambda: open_file())
file_choose_button.place(x = 20, y = 10)

#file window
filedir = Entry(my_window, text=" ", width = 100)#leave text empty
                                                 #as a place holder
filedir.place(x = 200, y = 10)

###############################################################################

L1 = Label(my_window, text = 'select level of significance')
L1.place(x = 20, y = 40)

L2 = Label(my_window, text = 'Sigma')
L2.place(x = 350, y = 40)

CB1 = ttk.Combobox(my_window, values = [1, 2, 3, 4, 5, 6])
CB1.place(x = 200, y = 40)

L3 = Label(my_window,text = 'select TR sample')
L3.place(x = 20, y = 70)

CB2 = ttk.Combobox(my_window, values = [1, 2, 3])
CB2.place(x = 200, y = 70)

###############################################################################

run_button = Button(my_window, text = 'run analysis',\
                            command = lambda: run())
run_button.place(x = 20, y = 100)

###############################################################################



def open_file():
    filename = askopenfilename()
    filedir.delete(0, "end")
    filedir.insert(0, filename)

def get_filename():
    filename = filedir.get()
    return(filename)

def get_sigma():
    sigma  = CB1.get()
    return(sigma)

def get_sample():
    sample = CB2.get()
    return(sample)

def run():
    filename = get_filename()
    sigma    = get_sigma()
    sample   = get_sample()
    
    sigma  = int(sigma)
    sample = int(sample)
    
    my_data_analysis(filename,sample,sigma)


my_window.mainloop()#runs an infinite loop