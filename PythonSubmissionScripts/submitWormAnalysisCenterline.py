#!/usr/bin/python
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
    import tkFileDialog
except ImportError:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog
import paramiko
import re
import math
import datetime

# need to make the master Tk window at the very beginning
master = Tk()


# store where we start from
startingDir = os.getcwd()
print((time.asctime()+'\n'))

# get ready for pickled variables 
pickle_path = (os.environ['HOME'] + "/platypusTemp/")
pickle_file = pickle_path + "pickles2.p"
if not os.path.exists(pickle_path):
	os.makedirs(pickle_path)
    
if not os.path.exists(pickle_file):
	storedUsername = { "username": "USER" }
	pickle.dump( storedUsername, open(pickle_file, "wb" ) )

# check to see if there is a username in the pickle file
prevUser = pickle.load( open( pickle_file, "rb" ) )
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
master.columnconfigure(7, pad=3)
master.columnconfigure(8, pad=3)
master.columnconfigure(9, pad=3)
master.columnconfigure(10, pad=3)

master.rowconfigure(0, pad=3)
master.rowconfigure(1, pad=3)

# user name
master.L1 = Label(master, text="User Name")
master.L1.grid(row=0, column=0, sticky=W+E)

master.e1 = Entry(master)
master.e1.insert(0, defaultName)
master.e1.grid(row=0, column=1, sticky=W+E)

# memory per job
master.L2 = Label(master, text="Memory per cell (mb)")
master.L2.grid(row=1, column=0, sticky=W+E)

master.e2 = Entry(master)
master.e2.insert(0,'20000')
master.e2.grid(row=1, column=1, sticky=W+E)

# time per cell
master.L3 = Label(master, text="time requested (min)")
master.L3.grid(row=2, column=0, sticky=W+E)

master.e3 = Entry(master)
master.e3.insert(0,'600')
master.e3.grid(row=2, column=1, sticky=W+E)

# path parent, normally tigress/LEIFER
master.L4 = Label(master, text="Parent Path")
master.L4.grid(row=3, column=0, sticky=W+E)

master.e4 = Entry(master)
master.e4.insert(0,'/tigress/LEIFER/PanNeuronal/')
master.e4.grid(row=3, column=1, sticky=W+E)

# data date
master.L5 = Label(master, text="Date of data")
master.L5.grid(row=4, column=0, sticky=W+E)

master.e5 = Entry(master)
master.e5.insert(0,defaultDate)
master.e5.grid(row=4, column=1, sticky=W+E)

# data folder
master.L6 = Label(master, text="DataFolderName")
master.L6.grid(row=5, column=0, sticky=W+E)

master.e6 = Entry(master)
master.e6.insert(0,defaultFolder)
master.e6.grid(row=5, column=1, sticky=W+E)



if os.name == 'posix':
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        
def submitScript(master=None):
    # username
    username = master.e1.get()
    print("Username: " + username)
    
    
    
    
    # memory requested
    memReq = master.e2.get()
    print("Memory requested: " + memReq + " mb per cell.")

    #time requested
    time = master.e3.get()
    
    server ='della.princeton.edu'
    print("Request: Use della")
    beginOfPath=master.e4.get()
    date=master.e5.get()
    brainFolder=master.e6.get()
    folderName=brainFolder
    isUpdateCode = 0 # no updating code or uploading data for now
    uploadCells = 0
    # which folder to process
    fullPath = os.path.join(beginOfPath,date)
    fullPath = os.path.join(fullPath, brainFolder)
    #picle dump
    prevUser['userName']=username
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
    
    
    
    print("full path is: " + fullPath)
    if isUpdateCode == 1:
        # see stackoverflow http://stackoverflow.com/questions/13851846/paramiko-sftpclient-setting-missing-host-key-policy
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy())
        try:
            client.connect('cetus.princeton.edu', 22, username)
        except paramiko.AuthenticationException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
        except SSHException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
            
        # untar code once it is on the server
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
            client.connect('cetus.princeton.edu', 22, username)
        except paramiko.AuthenticationException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
        except SSHException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
      
        commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
        commandList.insert(len(commandList)-1, "export GRID_HOME=$HOME")
        commandList.insert(len(commandList)-1, "cd scripts/shae-pythonSubmissionScripts")
        commandList.insert(len(commandList)-1, "echo environment variables")
        commandList.insert(len(commandList)-1, "echo $GRID_HOME")
        commandList.insert(len(commandList)-1, "python updateCodeFromGit.py")
        commands = "\n".join(commandList)
        #print(commands)
        stdin, stdout, stderr = client.exec_command(commands)
        print('Updating code from Git')
        returnedOutput = stdout.readlines()
        print(' '.join(returnedOutput))
        print('stdError: updating from Git')
        print(stderr.readlines())
        client.close()
        print('Done updating code.\n\n')
        
        
    if uploadCells == 1: 
        os.chdir(beginOfPath)
        
        # copy images over via sftp
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        try:
         client.connect(server, 22, username)
        except paramiko.AuthenticationException:
         password = master.e7.get()
         client.connect(server, 22, username, password)
        except SSHException:
         password = master.e7.get()
         client.connect(server, 22, username, password)
         
        client2 = client.open_sftp()
        client2.chdir('.') # change the current director to itself, not sure why this is necessary
        # move the file, hopefully to the right location
        try:
            client2.stat('data/' + folderName)
        except:
            client2.mkdir('data/' + folderName)

        print(os.path.join(folderName, fileName))
        client2.put(os.path.join(fullPath, fileName),("data/"+os.path.join(folderName, fileName)))
                            
        print('Done sending images.\n\n')

        
        #
        #commandList = ["pwd","pwd"]; # pwd at both ends, give the list something to add to the middle of
        #commandList.insert(len(commandList)-1, ("lcd "+startingDir));
        #commandList.insert(len(commandList)-1, "cd data");
        #commandList.insert(len(commandList)-1, "put tarredCells.tar.gz");
        #commands = "\n".join(commandList);
        #stdin, stdout, stderr = client.exec_command(commands)
        client2.close()
        client.close()
        print('Done transferring images.\n\n')

    # deal with folder names that have spaces
    matlabDirName =fullPath
    jobName = "a" + date
    # matlabDirName = "\\ ".join(matlabDirName)
    print(matlabDirName)
    
    # connect and submit job via qsub
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        client.connect(server, 22, username)
    except paramiko.AuthenticationException:
        password = master.e7.get()
        client.connect(server, 22, username, password)
    except SSHException:
        password = master.e7.get()
        client.connect(server, 22, username, password)
    
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    # set up the environment so that it matches an ssh login instead of the reduced paramiko one, hopefully this will help.
    commandList.insert(len(commandList)-1, "export SGE_ROOT=/sge")
    commandList.insert(len(commandList)-1, "export PATH=$PATH:/usr/local/emboss/bin:/sge/bin/lx-amd64:/usr/kerberos/bin")
    commandList.insert(len(commandList)-1, "export GRID_HOME=$HOME")
    commandList.insert(len(commandList)-1, "cd data")
    
    # make sure permissions on shell scripts are correct, they keep getting off for some reason
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormCellTracking.sh")
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormTrackCompiler.sh")
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormBotChecker.sh")
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormCenterlineFitting.sh")
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormCenterlineCompile.sh")
    commandList.insert(len(commandList)-1, "chmod 777 " + matlabDirName)

    # add somewhere for err and out files to go
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    commandList.insert(len(commandList)-1, "mkdir " + outputFilePath)
    
    #deal with folders with spaces in names

    matlabDirName = folderName
    jobName = folderName

   # matlabDirName = "\\ ".join(matlabDirName);
    print(folderName)

    print('Writing inputs line to text file')
    fileOutputName= fullPath + "/cmdOutput" +str(currentDate) +".txt"
    #fileOutputName= fullPath + "/cmdOutput_centerline" +str(currentDate) +".txt"
    #f=open(fileOutputName,'w')

    qString = " --time=" + time
    memString= " --mem=" + memReq
    userEmail="jnguyen@princeton.edu"
    jobName=folderName.split('_')[1]
    jobName='B'+jobName
#sbatch --mem=204800 --qos=1hr -n 4 -N 2-4 -D -o data/slurmtest -e data/slurmtest 
#--array=1-10:1 scripts/shae-pythonSubmissionScripts/runWormStraighten.sh 
#"20150616/BrainScanner20150616_100353-worm1"
    qsubCommand1 = ("sbatch" + memString + qString + " -D " + folderName
        + " -J "+ jobName
        + " --output=\"" + outputFilePath + "/CLjob-%J.out" + "\" --error=\"" 
        + outputFilePath + "/CLjob-%J.err" + "\""
        + " --array=1-32:1"
        + " --mail-type=end" + " --mail-user=" + userEmail
        + " ../scripts/shae-pythonSubmissionScripts/runWormCenterlineFitting.sh '" 
        + fullPath +"' ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch" + memString + qString + " -D " + folderName
        + " -J "+ jobName
        + " -d singleton"
        + " --output=\"" + outputFilePath + "/CLCompile-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/CLCompile-%J.err" + "\""
        + " --mail-type=end" + " --mail-user=" + userEmail
        + " ../scripts/shae-pythonSubmissionScripts/runWormCenterlineCompile.sh '" 
        + fullPath +"' ")
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    
    #f.write(qsubCommand1)
    #f.write('\n')

    commands = "\n".join(commandList)
    stdin, stdout, stderr = client.exec_command(commands)
    print('stdOutput: submitting job')
    returnedOutput = stdout.readlines()
    print(' '.join(returnedOutput))
    print('stdError: submitting job')
    print(stderr.readlines())
    client.close()
    print('Done submitting job.\n\n')
    print(matlabDirName)

    # close window at the end
    master.destroy()
        
                
def callback1(event=None,master=None):
    "Check for password and continue"
    print(master.e1.get())
    
    def enterPass1(event=None, master=None):
        "Request password to continue"
        if event == None:
            pass
        else:
            master = event.widget.master
        
        username = master.e1.get()
        print(username)
        # password, just the number of characters
        password = master.e7.get()
        
        passList = ["*"] * len(password)
        print("Entered password: " +   "".join(passList))
        
        submitScript(master)
    
    username = master.e1.get()
    if username == "noPass":
        isNeedsPassword = True
    else:
        # try to see if a password is needed
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        username = master.e1.get()
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
    
# bind enter key and button
master.b = Button(master, text="Enter", width=10, command=lambda:callback1(master=master))
master.b.grid(row=10,columnspan=2, sticky=W+E)
master.bind("<Return>", callback1b)

master.e1.focus_set()

mainloop()
