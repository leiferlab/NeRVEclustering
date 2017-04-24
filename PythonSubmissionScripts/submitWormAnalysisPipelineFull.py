#!/usr/bin/python



#code for submitting jobs onto Princeton Slurm cluster, mainly using Della now.
#Running this opens a GUI for selecting input files and parameters for analyzing
#imaging data after centerline and straightening steps. Neuron tracking algorithm
#is outlined in Nguyen et all (2016). This is the main part of the algorithm, with
#neuron registration vector encoding. Defaults work for almost all parameters now. 




import pickle
import os
import sys
import subprocess
import tarfile
import time
import glob
try:
    from Tkinter import *
    import tkFileDialog
except ImportError:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog

import paramiko
import datetime
import slurmInput as slurm
import socket 

# need to make the master Tk window at the very beginning
def make_gui():
    master = Tk()
    
    # store where we start from
    startingDir = os.getcwd()
    print((time.asctime()+'\n'))
    
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
    
    if 'frameNumber' in prevUser:
            defaultFrameNumber = prevUser['frameNumber']
    else:
            defaultFrameNumber = "1000"
            
    if 'refNumber' in prevUser:
            defaultRefNumber = prevUser['refNumber']
    else:
            defaultRefNumber = "100"

    if 'neuronNumber' in prevUser:
            defaultNeuronNumber = prevUser['neuronNumber']
    else:
            defaultNeuronNumber = "150"
            
    if 'checkNumber' in prevUser:
            defaultCheckNumber = prevUser['checkNumber']
    else:
            defaultCheckNumber = "500"
            
            
    # use the example of the calculator from http://zetcode.com/gui/tkinter/layout/ to help layout
    master.title("Options")
    ncol , nrow = 2 , 19
    for i in range(nrow):
        master.columnconfigure(i, pad=3)
    for i in range(ncol):
        master.rowconfigure(i, pad=3)
    
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
    
    # Number of frames
    L_nframe = Label(master, text="number of frame")
    L_nframe.grid(row=6, column=0, sticky=W+E)
    
    master.e['nframes'] = Entry(master)
    master.e['nframes'].insert(0,defaultFrameNumber)
    master.e['nframes'].grid(row=6, column=1, sticky=W+E)
    
    # Number of groups per frame
    L_ref = Label(master, text="Number of Reference")
    L_ref.grid(row=7, column=0, sticky=W+E)
    
    master.e['n_ref'] = Entry(master)
    master.e['n_ref'].insert(0,defaultRefNumber)
    master.e['n_ref'].grid(row=7, column=1, sticky=W+E)
    
    # Number of cells to check
    L_nneurons = Label(master, text="Number of Neurons")
    L_nneurons.grid(row=8, column=0, sticky=W+E)
    
    master.e['n_neurons'] = Entry(master)
    master.e['n_neurons'].insert(0,defaultNeuronNumber)
    master.e['n_neurons'].grid(row=8, column=1, sticky=W+E)
    
    # number fo frames to use from cross checking
    L_nchecks = Label(master, text="Number of checks")
    L_nchecks.grid(row=9, column=0, sticky=W+E)
    
    master.e['n_checks'] = Entry(master)
    master.e['n_checks'].insert(0,defaultCheckNumber)
    master.e['n_checks'].grid(row=9, column=1, sticky=W+E)


    #Run Straightening?
    L_straightflag = Label(master, text="Run Straightening")
    L_straightflag.grid(row=10, column=0, sticky=W+E)
    
    var1= IntVar()
    master.e['straight_flag'] = Checkbutton(master, text=None, variable=var1)
    master.e['straight_flag'].var = var1
    master.e['straight_flag'].grid(row=10, column=1, sticky=W+E)
    master.e['straight_flag'].var.set(1)
    
        
    # Run Tracker?
    L_trackflag = Label(master, text="Run Track")
    L_trackflag.grid(row=11, column=0, sticky=W+E)
    
    var2= IntVar()
    master.e['track_flag'] = Checkbutton(master, text=None, variable=var2)
    master.e['track_flag'].var = var2
    master.e['track_flag'].grid(row=11, column=1, sticky=W+E)
    master.e['track_flag'].var.set(1)
    
    
    # Run Checker?
    L_checkflag = Label(master, text="Run Check")
    L_checkflag.grid(row=13, column=0, sticky=W+E)
    
    var3= IntVar()
    master.e['check_flag'] = Checkbutton(master, text=None, variable=var3)
    master.e['check_flag'].var = var3
    master.e['check_flag'].grid(row=13, column=1, sticky=W+E)
    master.e['check_flag'].var.set(1)
    
    
    # Run cropper?
    L_cropflag = Label(master, text="Crop")
    L_cropflag.grid(row=15, column=0, sticky=W+E)
    
    var5= IntVar()
    master.e['crop_flag']= Checkbutton(master, text=None, variable=var5)
    master.e['crop_flag'].var = var5
    master.e['crop_flag'].grid(row=15, column=1, sticky=W+E)
    master.e['crop_flag'].var.set(1)
    
    L_email = Label(master, text="Crop")
    L_email.grid(row=16, column=0, sticky=W+E)
    
    var6= IntVar()
    master.e['mail_flag']= Checkbutton(master, text=None, variable=var5)
    master.e['mail_flag'].var = var6
    master.e['mail_flag'].grid(row=16, column=1, sticky=W+E)
    master.e['mail_flag'].var.set(1)
    
    return master

def get_nframes(username,fullPath):
    subDataFile = fullPath + '/submissionParameters.txt'
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    client.connect(server,22,username)
    client2 = client.open_sftp()
    remote_file=client2.open(subDataFile)
    sub_data={}
    for lines in remote_file:
        lines=lines.split(" ")
        print(lines)
        sub_data[lines[0]]=lines[1]
        
    if sub_data.get('NFrames') != None:
        print('Getting number of frames from Remote')
        NFrames=sub_data['NFrames']
        return NFrames
    else: 
        raise NameError("Cannot connect to get number of Frames")
    

        
        
def submitScript(master=None):
    server ='della.princeton.edu'
    # username
    username = master.e['user_name'].get()
    print("Username: " + username)    
    
    #construct path
    beginOfPath=master.e['parent_path'].get()
    date=master.e['date'].get()
    folderName=master.e['folder_name'].get()
    
    # which folder to process, must add paths linux style
    fullPath = beginOfPath + "/" + date + "/" + folderName
    
    #read in N frames using sftp, or else, get from GUI
    try:
        NFrames=get_nframes(username,fullPath)
    except:
        #at some point I'll learn how to do proper exceptions
        print('Getting number of frames from GUI')
        NFrames=master.e['nframes'].get()
            
    print("Number of Files is" + NFrames)
    totalRuns=str(int(NFrames))
    print("Total runs is " + totalRuns)

    #number of ref frames
    nRef=master.e['n_ref'].get()
    
    # groups  for bot checker
    nNeurons=master.e['n_neurons'].get()
    print("Number of Neruons is " + nNeurons)
    
    nCheck=master.e['n_checks'].get()
    
    #save defaults using pickle dump
    pickle_path = (os.path.expanduser('~') + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    prevUser=slurm.pickle_load()
    prevUser['username']=username
    prevUser['frameNumber']=NFrames
    prevUser['date'] = date
    prevUser['folderName']=folderName
    prevUser['refNumber']=nRef
    prevUser['neuronNumber']=nNeurons
    prevUser['checkNumber']=nCheck
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
    
    totalRuns=int(totalRuns)
    totalTime=270
    totalTime=str(totalTime)
  
    #Which parts of code to run
    straightFlag    = master.e['straight_flag'].var.get()
    trackFlag       = master.e['track_flag'].var.get()
    checkFlag       = master.e['check_flag'].var.get()
    cropFlag        = master.e['crop_flag'].var.get()
    emailFlag       = master.e['email_flag'].var.get()
    
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
        password = master.e['nframes'].get()
        client.connect(server, 22, username, password)
    except SSHException:
        password = master.e['nframes'].get()
        client.connect(server, 22, username, password)
    
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    # set up the environment so that it matches an ssh login instead of the reduced paramiko one, hopefully this will help.
        
    # add somewhere for err and out files to go
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    commandList.insert(len(commandList)-1, "mkdir " + outputFilePath)
    
    #deal with folders with spaces in names

   # matlabDirName = "\\ ".join(matlabDirName);
    print('Writing inputs line to text file')
    userEmail=username+"@princeton.edu"
    
    #submit path setup bash commands
    commandList=slurm.path_setup(commandList)
    #make job name from foldername
    if straightFlag:
        commandList=slurm.straighten_input(commandList,fullPath,totalRuns)
        
    if trackFlag:
        commandList=slurm.track_input(commandList,fullPath,totalRuns,nRef)
        
    if checkFlag:
        commandList=slurm.check_input(commandList,fullPath,totalRuns,nCheck,nNeurons)
        
    if cropFlag:
        commandList=slurm.crop_input(commandList,fullPath,emailFlag)

    commands = "\n".join(commandList)
    print(commands)
    
    slurm.make_ouputfolder(client,fullPath)
    #write commands to text file via paramiko
    slurm.write_input(commandList,client,fullPath)
    
    stdin, stdout, stderr = client.exec_command(commands)
    print('stdOutput: submitting job')
    returnedOutput = stdout.readlines()
    print(' '.join(returnedOutput))
    print('stdError: submitting job')
    print(stderr.readlines())
    client.close()
    print('Done submitting job.\n\n')

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
    password = master.e['nframes'].get()

    passList = ["*"] * len(password)
    print("Entered password: " +   "".join(passList))

    submitScript(master)
                
def callback1(event=None,master=None):
    "Check for password and continue"
    print(master.e['user_name'].get())
    
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
        username = master.e['user_name'].get()
        try:
            client.connect('della.princeton.edu', 22, username,pkey=key)
            isNeedsPassword = False
        except paramiko.AuthenticationException:
            isNeedsPassword = True
            
    if isNeedsPassword:
        # use the same window as before, just add an additional password field
        # password
        master.L15 = Label(master, text="Password")
        master.L15.grid(row=14, column=0, sticky=W+E)
        
        master.e15 = Entry(master,show="*")
        master.e15.insert(0, "Password")
        master.e15.grid(row=14, column=1, sticky=W+E)
        
        master.b = Button(master, text="Enter", width=10, command=lambda:enterPass1(master=master))
        master.b.grid(row=16,columnspan=2, sticky=W+E)
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
        This is the submission script for running analysis on whole brain imaging.
        This code runs straightening, tracking, and cross validation of points. It also 
        extracts signals and produces heatmaps. This can be run after the wormAnalysisPreview
        gui gives all green lights. Check the readme for instructions for della access. 
        
        
        For a quick test, run this code as follows:
        User Name: <your username>
        Parent Path:/tigress/LEIFER/PanNeuronal
        Date of Data: testing_sets
        Data Folder Name: BrainScanner20161031_111303
        Number of Frames : 1000
        Number of References: 10
        Number of Neurons: 150
        Number of Checks: 100
        <click all check boxes>
        
        ''')
    master=make_gui()
    master.b = Button(master, text="Enter", width=10, command=lambda:callback1(master=master))
    master.b.grid(row=17,columnspan=2, sticky=W+E)
    master.bind("<Return>", callback1b)
    
    if  socket.gethostname()=='tigressdata.princeton.edu':
        master.b2 = Button(master, text='Select Folder', width=10, command=lambda:SelectFolder(master=master))
        master.b2.grid(row=18,columnspan=2, sticky=W+E)
        
    master.e['user_name'].focus_set()
    
    if os.name == 'posix':
        os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        
    mainloop()
