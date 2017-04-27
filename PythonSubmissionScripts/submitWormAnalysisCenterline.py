#!/usr/bin/python
# 
#code for submitting jobs onto Princeton Slurm cluster, mainly using Della now.
#Running this opens a GUI for selecting input files and parameters for analyzing
#imaging data after centerline and straightening steps. Neuron tracking algorithm
#is outlined in Nguyen et all (2016). This is the Centerline portion of the code. 
#The folder requires an alignment. mat file
# 
import pickle
import os
import time
try:
    from Tkinter import *
    import tkFileDialog
except ImportError:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog

import paramiko
import slurmInput as slurm
import socket


def make_gui():
    # need to make the master Tk window at the very beginning
    master = Tk()
    
    # store where we start from
    startingDir = os.getcwd()
    print((time.asctime()+'\n'))
    
    # get ready for pickled variables 
    prevUser=slurm.pickle_load()
    
    if 'username' in prevUser:
    	defaultName = prevUser['username']
    else:
    	defaultName = "USER"
        
    if 'date' in prevUser:
    	defaultDate = prevUser['date']
    else:
    	defaultDate = "20160101"
            
    if 'folderName' in prevUser:
    	defaultFolder = prevUser['folderName']
    else:
    	defaultFolder = "Brain_Scanner2016"
    
             
    # use the example of the calculator from http://zetcode.com/gui/tkinter/layout/ to help layout
    master.title("Options")
    
    master.columnconfigure(0, pad=3)
    master.columnconfigure(1, pad=3)
    master.columnconfigure(2, pad=3)
    master.columnconfigure(3, pad=3)
    master.columnconfigure(4, pad=3)
    master.columnconfigure(5, pad=3)
    master.columnconfigure(6, pad=3)
    master.columnconfigure(9, pad=3)
    master.columnconfigure(10, pad=3)
    master.columnconfigure(11, pad=3)
    
    master.rowconfigure(0, pad=3)
    master.rowconfigure(1, pad=3)
    
    master.e=dict()
    
    # user name
    L_username = Label(master, text="User Name")
    L_username.grid(row=0, column=0, sticky=W+E)
    
    master.e['user_name'] = Entry(master)
    master.e['user_name'].insert(0, defaultName)
    master.e['user_name'].grid(row=0, column=1, sticky=W+E)
    
    
    # path parent, normally tigress/LEIFER
    L_path = Label(master, text="Parent Path")
    L_path.grid(row=3, column=0, sticky=W+E)
    
    master.e['parent_path'] = Entry(master)
    master.e['parent_path'].insert(0,'/tigress/LEIFER/PanNeuronal')
    master.e['parent_path'].grid(row=3, column=1, sticky=W+E)
    
    # data date
    L_date = Label(master, text="Date of data")
    L_date.grid(row=4, column=0, sticky=W+E)
    
    master.e['date'] = Entry(master)
    master.e['date'].insert(0,defaultDate)
    master.e['date'].grid(row=4, column=1, sticky=W+E)
    
    # data folder
    L_foldername = Label(master, text="DataFolderName")
    L_foldername.grid(row=5, column=0, sticky=W+E)
    
    master.e['folder_name'] = Entry(master)
    master.e['folder_name'].insert(0,defaultFolder)
    master.e['folder_name'].grid(row=5, column=1, sticky=W+E)
    
    L_email = Label(master, text="Email")
    L_email.grid(row=6, column=0, sticky=W+E)
    
    var1= IntVar()
    master.e['mail_flag']= Checkbutton(master, text=None, variable=var1)
    master.e['mail_flag'].var = var1
    master.e['mail_flag'].grid(row=6, column=1, sticky=W+E)
    master.e['mail_flag'].var.set(1)
    
    return master
    

def submitScript(master=None):
    # username
    username = master.e['user_name'].get()
    print("Username: " + username)
    
    
    server ='della.princeton.edu'
    print("Request: Use della")
    beginOfPath=master.e['parent_path'].get()
    date=master.e['date'].get()
    folderName=master.e['folder_name'].get()
    print("dateFolder is " + date)
    print("beginOfPath is "+ beginOfPath)
    print(folderName)
    # which folder to process, must add paths linux style
    fullPath = beginOfPath + "/" + date
    fullPath = fullPath + "/" + folderName
    
    #picle dump
    #save defaults using pickle dump
    pickle_path = (os.environ['HOME'] + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    prevUser=slurm.pickle_load()
    prevUser['username']=username
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
    print("full path is: " + fullPath)
    # deal with folder names that have spaces
    
    if socket.gethostname()=='tigressdata.princeton.edu':
        keypath='/tigress/LEIFER/.ssh/id_rsa'
        key = paramiko.RSAKey.from_private_key_file(keypath)
    else:
        key=None
    # connect and submit job via sbatch
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        client.connect(server, 22, username,pkey=key)
    except paramiko.AuthenticationException:
        password = master.e7.get()
        client.connect(server, 22, username, password)
    except SSHException:
        password = master.e7.get()
        client.connect(server, 22, username, password)
    
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    # set up the environment so that it matches an ssh login instead of the reduced paramiko one, hopefully this will help.
    
    # add somewhere for err and out files to go
    outputFilePath=slurm.make_output_path(fullPath)
    stdin,stdout,stderr = client.exec_command("mkdir -p " + outputFilePath)
    
    
    print('Writing inputs line to text file')

    commandList=slurm.path_setup(commandList)
    commandList=slurm.centerline_input(commandList,fullPath)
    
    commands = "\n".join(commandList)
    stdin, stdout, stderr = client.exec_command(commands)
    slurm.write_input(commandList,client,fullPath)

    print('stdOutput: submitting job')
    
    returnedOutput = stdout.readlines()
    print(' '.join(returnedOutput))
    print('stdError: submitting job')
    print(stderr.readlines())
    client.close()
    print('Done submitting job.\n\n')

    print('''
        Output files will be saved in 
        '''
        + fullPath
        + '''
        behaviorAnalysis Folder with a centerline.mat inside
        ''')
    # close window at the end
    master.destroy()
        

def enterPass1(event=None, master=None):
    "Request password to continue"
    if event == None:
        pass
    else:
        master = event.widget.master
    
    username = master.e['user_name'].get()
    print(username)
    # password, just the number of characters
    password = master.e7.get()
    
    passList = ["*"] * len(password)
    print("Entered password: " +   "".join(passList))
    
    submitScript(master)
    

def callback1(event=None,master=None):
    "Check for password and continue"
    print(master.e['user_name'] .get())
    
    
    username = master.e['user_name'].get()
    if username == "noPass":
        isNeedsPassword = True
    else:
        if socket.gethostname()=='tigressdata.princeton.edu':
            keypath='/tigress/LEIFER/.ssh/id_rsa'
            key = paramiko.RSAKey.from_private_key_file(keypath)
        else:
            key=None
        # try to see if a password is needed
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        username = master.e['user_name'] .get()
        try:
            client.connect('della.princeton.edu', 22, username)
            isNeedsPassword = False
        except paramiko.AuthenticationException:
            isNeedsPassword = True
            
    if isNeedsPassword:
        # use the same window as before, just add an additional password field
        # password
        master.L6 = Label(master, text="Password")
        master.L6.grid(row=6, column=0, sticky=W+E)
        
        master.e7 = Entry(master,show="*")
        master.e7.insert(0, "Password")
        master.e7.grid(row=7, column=1, sticky=W+E)
        
        master.b = Button(master, text="Enter", width=10, command=lambda:enterPass1(master=master))
        master.b.grid(row=6,columnspan=2, sticky=W+E)
        master.bind("<Return>", enterPass1) 
    else:
        print("No password needed")
        submitScript(master)
        
        
def callback1b(event=None):
    "Grabs the master of whatever widget exists in event and executes the appropriate callback"
    master = event.widget.master
    #master.withdraw()
    callback1(master=master)
    
def SelectFolder(master=None):
   
    folder=tkFileDialog.askdirectory(mustexist=False , initialdir= '/tigress/LEIFER/PanNeuronal/')
    if folder:
        path,folderName=os.path.split(folder)
        path,date=os.path.split(path)
        master.e['parent_path'].delete(0,END)
        master.e['parent_path'].insert(0,path)
        master.e['date'].delete(0,END)
        master.e['date'].insert(0,date)
        master.e['folder_name'].delete(0,END)
        master.e['folder_name'].insert(0,folderName)
        print folder
    
if __name__ == '__main__':
# bind enter key and button
    print('''
        This is the submission script for running time alignment and detcting flashes in videos. 
        The lowMag folders must be inside the corresponding BrainScanner folder on della. 
        The videos must each contain at least 2 flashes.
        
        
        For a quick test, run this code as follows:
        User Name: <your username>
        Parent Path:/tigress/LEIFER/PanNeuronal
        Date of Data: testing_sets
        Data Folder Name: BrainScanner20161031_111303
        
        ''')
    master=make_gui()
    
    # bind enter key and button
    master.b = Button(master, text="Enter", width=10, command=lambda:callback1(master=master))
    master.b.grid(row=10,columnspan=2, sticky=W+E)
    master.bind("<Return>", callback1b)


    if  socket.gethostname()=='tigressdata.princeton.edu':
        master.b2 = Button(master, text='Select Folder', width=10, command=lambda:SelectFolder(master=master))
        master.b2.grid(row=11,columnspan=2, sticky=W+E)
        
    master.e['user_name'].focus_set()
    
    mainloop()
