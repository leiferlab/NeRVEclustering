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

if 'frameNumber' in prevUser:
	defaultFrameNumber = prevUser['frameNumber']
else:
	defaultFrameNumber = "1000"
    

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
master.columnconfigure(11, pad=3)
master.columnconfigure(12, pad=3)
master.columnconfigure(13, pad=3)
master.columnconfigure(14, pad=3)
master.columnconfigure(15, pad=3)

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
master.e2.insert(0,'8000')
master.e2.grid(row=1, column=1, sticky=W+E)

# time per cell
master.L3 = Label(master, text="time requested (min)")
master.L3.grid(row=2, column=0, sticky=W+E)
master.e3 = Entry(master)
master.e3.insert(0,'240')
master.e3.grid(row=2, column=1, sticky=W+E)

# path parent, normally /tigress/LEIFER
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


#Slurm or SGE? Default use slurm
master.L7 = Label(master, text="Which Cluster?")
master.L7.grid(row=6, column=0, sticky=W+E)

var3 = StringVar()
master.CB7 = OptionMenu(master, var3, "Della","Genomics","SGE")
master.CB7.var3 = var3
master.CB7.grid(row=6, column=1, sticky=W+E)
master.CB7.var3.set("Della")

master.L8 = Label(master, text="Number of Volumes")
master.L8.grid(row=7, column=0, sticky=W+E)
master.e8 = Entry(master)
master.e8.insert(0,defaultFrameNumber)
master.e8.grid(row=7, column=1, sticky=W+E)

master.L9 = Label(master, text="Step size")
master.L9.grid(row=8, column=0, sticky=W+E)
master.e9 = Entry(master)
master.e9.insert(0,'30')
master.e9.grid(row=8, column=1, sticky=W+E)

# Run Track compiler?
master.L10 = Label(master, text="Straighten Start")
master.L10.grid(row=9, column=0, sticky=W+E)

var0= IntVar()
master.CB10 = Checkbutton(master, text=None, variable=var0)
master.CB10.var = var0
master.CB10.grid(row=9, column=1, sticky=W+E)
master.CB10.var.set(1)


# Run Tracker?
master.L11 = Label(master, text="Run Straighten")
master.L11.grid(row=10, column=0, sticky=W+E)

var1= IntVar()
master.CB11 = Checkbutton(master, text=None, variable=var1)
master.CB11.var = var1
master.CB11.grid(row=10, column=1, sticky=W+E)
master.CB11.var.set(1)

# Run Track compiler?
master.L12 = Label(master, text="Straighten Compile")
master.L12.grid(row=11, column=0, sticky=W+E)

var2= IntVar()
master.CB12 = Checkbutton(master, text=None, variable=var2)
master.CB12.var = var2
master.CB12.grid(row=11, column=1, sticky=W+E)
master.CB12.var.set(1)

if os.name == 'posix':
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        
def submitScript(master=None):
    # username
    username = master.e1.get()
    print("Username: " + username)
        
    # memory requested
    memReq = master.e2.get()
    print("Memory requested: " + memReq + " gb per cell.")

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
    totalRuns=master.e8.get()
    stepSize=master.e9.get()
    cluster=master.CB7.var3.get()
    
    # which folder to process, must add paths linux style
    fullPath = beginOfPath + "/" + date
    fullPath = fullPath + "/" + brainFolder
    
    straightenFlag       = master.CB11.var.get()
    straightenCompileFlag= master.CB12.var.get()
    straightenStartFlag  = master.CB10.var.get()
    subDataFile = fullPath + '/submissionParameters.txt'
    try:
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
    
        # grouping size for step 1
        if sub_data.get('NFrames') != None:
            print('Getting number of frames from Remote')
            totalRuns=sub_data['NFrames']
        else:
            print('Getting number of frames from GUI')
            totalRuns=master.e8.get()
    except:
        #at some point I'll learn how to do proper exceptions
            print('Getting number of frames from GUI')
            totalRuns=master.e8.get()
    print("totalsRuns is "+ totalRuns)
    
    #save inputs
    prevUser['username']=username
    prevUser['frameNumber']=totalRuns
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
    
    if cluster in ["Genomics","Della"]:
       slurmFlag = 1
       totalRuns=float(totalRuns)
       nRuns=math.ceil(totalRuns/1000)
       # can't run more than 1000 jobs on slurm at once  =(
       totalTime=str(time)
       totalRuns=str(int(totalRuns))
       stepSize=str(stepSize)
       totalTime=str(int(stepSize)*5)
       if cluster=="Genomics":
           server ='gen-comp1.princeton.edu'
           dellaFlag=0
       elif cluster=="Della":
           server ='della.princeton.edu'
           dellaFlag = 1
       
    elif cluster=="Della":
       slurmFlag = 0
       server='cetus.princeton.edu'
       dellaFlag = 0

    print("Request: Use "+ cluster)

    isUpdateCode = 0
    
    if isUpdateCode == 1:
        # see stackoverflow http://stackoverflow.com/questions/13851846/paramiko-sftpclient-setting-missing-host-key-policy
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
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormStraighten.sh")
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormCompilePointStats.sh")
    commandList.insert(len(commandList)-1, "chmod 777 " + matlabDirName)

    # add somewhere for err and out files to go
    commandList.insert(len(commandList)-1, "mkdir " + folderName + "/outputFiles")
    # add somewhere for err and out files to go
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    commandList.insert(len(commandList)-1, "mkdir " + outputFilePath)
    
    #deal with folders with spaces in names

    matlabDirName = folderName
    jobName = folderName
    qString0 = "--time=50"
    qString1 = "--time=" + str(totalTime)
    qString2 = "--time=300" 
    qString3 = "--time=300"        
#sbatch --mem=204800 --qos=1hr -n 4 -N 2-4 -D -o data/slurmtest -e data/slurmtest 
#--array=1-10:1 scripts/shae-pythonSubmissionScripts/runWormStraighten.sh 
#"20150616/BrainScanner20150616_100353-worm1"
# 
# 

#    	qsubCommand1 = ("sbatch --mem=8000 " + qString2 + " -D " + folderName
#    		+ " -J "+ folderName + " -d singleton"
#    		+ " --output=\"" + "outputFiles/slurmJobStriaghtenInitial.out" + "\" --error=\"" 
#    		+ "outputFiles/slurmJobStriaghtenInitial.err \"" + " " 
#    		+ " ../scripts/shae-pythonSubmissionScripts/runWormStriaghtenInitial.sh '" 
#    		+ fileName +"'")
    if straightenStartFlag:
        input0 = "clusterStraightenStart('"+ fullPath + "')"
        qsubCommand0 = ("sbatch --mem=" + memReq +" "+ qString0 + " -D " + folderName
            + " -J "+ folderName
            + " --output=\"" + outputFilePath + "/CLstart-%J.out" + "\""
            + " --error=\"" + outputFilePath + "/CLstart-%J.err" + "\""
            + " ../scripts/shae-pythonSubmissionScripts/runMatlabInput.sh \"" 
            + input0 +"\"")
        qsubCommand0 = "q0=$("+ qsubCommand0 + ")"
        commandList.insert(len(commandList)-1, qsubCommand0)
        commandList.insert(len(commandList)-1, "echo $q0")
        print(qsubCommand0)
        dependencyString=" --dependency=afterok:${q0##* }"
    else:
        dependencyString=""
        
            
            
    if straightenFlag:
    	for i in range(0,int(nRuns)):
    		offset=str(i*1000)
    		if i>=(nRuns-1):
    			currentLimit=str(int(totalRuns)%1000)
    		else:
    			currentLimit="1000"

    		qsubCommand1 = ("sbatch --mem=" + memReq +" "+ qString1 + " -D " + folderName
    	        + " -J "+ folderName 
                + dependencyString
    	        + " --output=\"" + outputFilePath + "/slurmJob-%J.out" + "\" --error=\"" 
    	        + outputFilePath + "/slurmJob-%J.err" + "\""
    	        + " --array=1-" + currentLimit + ":" + str(stepSize) 
    	        + " ../scripts/shae-pythonSubmissionScripts/runWormStraighten.sh '" 
    	        + fullPath +"' "  + str(stepSize) +" " + offset)
    		commandList.insert(len(commandList)-1, qsubCommand1)
    		print(qsubCommand1)
            
    if straightenCompileFlag:
        qsubCommand2 = ("sbatch --mem=2000 " + qString2 + " -D " + folderName
            + " -J "+ folderName + " -d singleton"
            + " --output=\"" + outputFilePath + "/pscompile-%J.out" + "\" "
            + "--error=\"" + outputFilePath + "/pscompile-%J.err" + "\""
            + " ../scripts/shae-pythonSubmissionScripts/runWormCompilePointStats.sh '" 
            + fullPath +"'")
        print(qsubCommand2)
        commandList.insert(len(commandList)-1, qsubCommand2)


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
        client.set_missing_host_key_policy(paramiko.WarningPolicy())
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
master.b = Button(master, text="Enter", width=15, command=lambda:callback1(master=master))
master.b.grid(row=15,columnspan=2, sticky=W+E)
master.bind("<Return>", callback1b)

master.e1.focus_set()

mainloop()
