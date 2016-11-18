#code for submitting jobs onto Princeton Slurm cluster, mainly using Della now.
#Running this opens a GUI for selecting input files and parameters for analyzing
#imaging data after centerline and straightening steps. Neuron tracking algorithm
#is outlined in Nguyen et all (2016). This is the main part of the algorithm, with
#neuron registration vector encoding. Defaults work for almost all parameters now. 




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
master.columnconfigure(16, pad=3)

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
master.e2.insert(0,'2000')
master.e2.grid(row=1, column=1, sticky=W+E)

# time per cell
master.L3 = Label(master, text="time requested (min)")
master.L3.grid(row=2, column=0, sticky=W+E)

master.e3 = Entry(master)
master.e3.insert(0,'240')
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

# Number of frames
master.L7 = Label(master, text="number of frame")
master.L7.grid(row=6, column=0, sticky=W+E)

master.e7 = Entry(master)
master.e7.insert(0,defaultFrameNumber)
master.e7.grid(row=6, column=1, sticky=W+E)

# Number of groups per frame
master.L8 = Label(master, text="groups per frame")
master.L8.grid(row=7, column=0, sticky=W+E)

master.e8 = Entry(master)
master.e8.insert(0,'2')
master.e8.grid(row=7, column=1, sticky=W+E)

# Number of cells to check
master.L9 = Label(master, text="Number of Neurons")
master.L9.grid(row=8, column=0, sticky=W+E)

master.e9 = Entry(master)
master.e9.insert(0,'150')
master.e9.grid(row=8, column=1, sticky=W+E)

# Number of groups per cell
master.L10 = Label(master, text="Number of Groups per cell")
master.L10.grid(row=9, column=0, sticky=W+E)

master.e10 = Entry(master)
master.e10.insert(0,'1')
master.e10.grid(row=9, column=1, sticky=W+E)

# Run Tracker?
master.L11 = Label(master, text="Run Track")
master.L11.grid(row=10, column=0, sticky=W+E)

var1= IntVar()
master.CB11 = Checkbutton(master, text=None, variable=var1)
master.CB11.var = var1
master.CB11.grid(row=10, column=1, sticky=W+E)
master.CB11.var.set(0)

# Run Track compiler?
master.L12 = Label(master, text="Compile Track")
master.L12.grid(row=11, column=0, sticky=W+E)

var2= IntVar()
master.CB12 = Checkbutton(master, text=None, variable=var2)
master.CB12.var = var2
master.CB12.grid(row=11, column=1, sticky=W+E)
master.CB12.var.set(0)

# Run Checker?
master.L13 = Label(master, text="Run Check")
master.L13.grid(row=12, column=0, sticky=W+E)

var3= IntVar()
master.CB13 = Checkbutton(master, text=None, variable=var3)
master.CB13.var = var3
master.CB13.grid(row=12, column=1, sticky=W+E)
master.CB13.var.set(0)

# Run Check Compiler?
master.L14 = Label(master, text="Compile Check")
master.L14.grid(row=13, column=0, sticky=W+E)

var4= IntVar()
master.CB14= Checkbutton(master, text=None, variable=var4)
master.CB14.var = var4
master.CB14.grid(row=13, column=1, sticky=W+E)
master.CB14.var.set(0)

# Run cropper?
master.L15 = Label(master, text="Crop")
master.L15.grid(row=14, column=0, sticky=W+E)

var5= IntVar()
master.CB15= Checkbutton(master, text=None, variable=var5)
master.CB15.var = var5
master.CB15.grid(row=14, column=1, sticky=W+E)
master.CB15.var.set(0)


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
    beginOfPath=master.e4.get()
    date=master.e5.get()
    brainFolder=master.e6.get()
    folderName=brainFolder
    # which folder to process, must add paths linux style
    fullPath = beginOfPath + "/" + date
    fullPath = fullPath + "/" + brainFolder
    
    #read in N frames using sftp, or else, get from GUI
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
            NFrames=sub_data['NFrames']
        else:
            print('Getting number of frames from GUI')
            NFrames=master.e7.get()
    except:
        #at some point I'll learn how to do proper exceptions
            print('Getting number of frames from GUI')
            NFrames=master.e7.get()
            
    print("Number of Files is" + NFrames)
    nGroups=master.e8.get()
    print("Number of Groups is " + nGroups)
    totalRuns=str(int(NFrames)*int(nGroups))
    print("Total runs is " + totalRuns)

    # groups  for bot checker
    neuronRuns=master.e9.get()
    print("Number of Neruons is " + neuronRuns)
    groupsPerCell=master.e10.get()
    print("Number of groups per neuron is " + groupsPerCell)
    totalNeuronRuns=str(int(neuronRuns)*int(groupsPerCell))

    server ='della.princeton.edu'
    print("Request: Use della")

#    pickle.dump({"username": username, 
#        "frameNumber": NFrames,
#        "date": date, 
#        "folderName": folderName},
#         open(pickle_file, "wb" ) )
    prevUser['username']=username
    prevUser['frameNumber']=NFrames
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
    
    totalRuns=float(totalRuns)
    nRuns=math.ceil(totalRuns/1000)
       # can't run more than 1000 jobs on slurm at once  =(
    stepSize=max([1,int(math.ceil(totalRuns/800))])
    totalRuns=int(math.ceil(totalRuns/stepSize))
    totalTime=stepSize*90
    totalTime=str(totalTime)
    totalRuns=str(totalRuns)
  
    #Which parts of code to run
    trackFlag       = master.CB11.var.get()
    trackCompileFlag= master.CB12.var.get()
    checkFlag       = master.CB13.var.get()
    checkCompileFlag= master.CB14.var.get()
    cropFlag        = master.CB15.var.get()
    
    print("Request: Use "+ server)

    isUpdateCode = 0
    uploadCells=0
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
    fileName = 'PointsStats.mat'
    matlabDirName = os.path.join(fullPath, fileName)
    jobName = "_".join(matlabDirName)
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
    commandList.insert(len(commandList)-1, "chmod 777 ../scripts/shae-pythonSubmissionScripts/runWormBotCheckCompiler.sh")
    commandList.insert(len(commandList)-1, "chmod 777 " + matlabDirName)

    # add somewhere for err and out files to go
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    commandList.insert(len(commandList)-1, "mkdir " + outputFilePath)
    
    #deal with folders with spaces in names
    fileName2 = 'PointsStats2.mat'
    matlabDirName2 = fullPath + "/" + fileName2
   # matlabDirName = "\\ ".join(matlabDirName);
    print('Writing inputs line to text file')
    currentDate=datetime.date.today()
    slurmFlag=1
    dellaFlag=1;
    userEmail=username+"@princeton.edu"
    
    #make job name from foldername
    jobName=folderName.split('_')[1]
    jobName='B'+jobName
    if slurmFlag ==1:
        #genomics takes in --qos while della only wants estimated time and assigns queues itself
        if dellaFlag==0:
            qString1 = "--qos=1day"
            qString2 = "--qos=1day"
            qString3 = "--qos=1day"
        elif dellaFlag==1:
            qString1 = "--time=" + str(80*stepSize)
            qString2 = "--time=180" 
            qString3 = "--time=600"        

        if trackFlag:
            for i in range(0,int(nRuns)):
                #number of jobs capped at 1000, work around by submitting many batches with length 1000. 
                offset=str(i*1000)
                if i>=(nRuns-1):
                    currentLimit=str(int(totalRuns)%1000)
                else:
                    currentLimit="1000"
    
                qsubCommand1 = ("sbatch --mem=4000 " + qString1 + " -D " + folderName
                    + " -J "+ jobName
                    + " --output=\"" + outputFilePath + "/trackJob-%J.out" + "\"" 
                    + " --error=\"" + outputFilePath + "/trackJob-%J.err" + "\""
                    + " --array=1-" + currentLimit + ":" + str(stepSize)
                    + " ../scripts/shae-pythonSubmissionScripts/runWormCellTracking.sh '" 
                    + matlabDirName +"' " + nGroups +" " + offset +" " + str(stepSize))
                commandList.insert(len(commandList)-1, qsubCommand1)
                print(qsubCommand1)            

        if trackCompileFlag:
            qsubCommand2 = ("sbatch --mem=100000 " + qString2 + " -D " + folderName
                + " -J "+ jobName + " -d singleton"
                + " --output=\"" + outputFilePath + "/trackCompilerJob-%J.out" + "\""
                + " --error=\"" + outputFilePath + "/trackCompilerJob-%J.err" + "\""
                + " --mail-type=end" + " --mail-user=" + userEmail
                + " ../scripts/shae-pythonSubmissionScripts/runWormTrackCompiler.sh '" 
                + matlabDirName +"' '" + matlabDirName2 + "'")
            commandList.insert(len(commandList)-1, qsubCommand2)
            print(qsubCommand2)
            
        if checkFlag:
            qsubCommand3 = ("sbatch --mem=8000 " + qString3 + " -D " + folderName
                + " -J "+ jobName + " -d singleton"
                + " --output=\"" + outputFilePath + "/botCheckJob-%J.out" + "\""
                + " --error=\"" + outputFilePath + "/botCheckJob-%J.err" + "\""
                + " --array=1-"+ totalNeuronRuns + ":" + groupsPerCell
                + " ../scripts/shae-pythonSubmissionScripts/runWormBotChecker.sh '" 
                + matlabDirName2 +"' " + groupsPerCell + " 1")
            commandList.insert(len(commandList)-1, qsubCommand3)
            print(qsubCommand3)
                
        if checkCompileFlag:
            qsubCommand4 = ("sbatch --mem=100000 " + qString3 + " -D " + folderName
                + " -J "+ jobName + " -d singleton"
                + " --output=\"" + outputFilePath + "/botCompileJob-%J.out" + "\""
                + " --error=\"" + outputFilePath + "/botCompileJob-%J.err" + "\""
                + " --mail-type=end" + " --mail-user=" + userEmail
                + " ../scripts/shae-pythonSubmissionScripts/runWormBotCheckCompiler.sh  '" 
                + fullPath +"' ")
            commandList.insert(len(commandList)-1, qsubCommand4)
            print(qsubCommand4)
            
        if cropFlag:
            input1= "fiducialCropper3('"+ fullPath +"')"
            print(input1)
            qsubCommand5 = ("sbatch --mem=16000 " + qString2 + " -D " + folderName
                + " -J "+ jobName + " -d singleton"
                + " --output=\"" + outputFilePath + "/crop-%J.out"+ "\"" 
                + " --error=\"" + outputFilePath + "/crop-%J.err" + "\""
                + " --mail-type=end" + " --mail-user=" + userEmail
                + " ../scripts/shae-pythonSubmissionScripts/runMatlabInput.sh \"" 
                + input1 +"\" ")
            commandList.insert(len(commandList)-1, qsubCommand5)
            print(qsubCommand5)
    else:
    #SGE has a problem with memory, asking for 20g of mem doesnt actually do anything, 
        qsubCommand1 = ("qsub" + " -l " + queueUsed + " -l gb=" + memReq + " -l h_data=" 
                + memReq + "g -cwd" + " -o \"" + folderName + "\" -e \"" 
                + folderName + "\" -t 1-" + totalRuns + ":" + str(1) 
                + " -N \"" + jobName + "\" ../scripts/shae-pythonSubmissionScripts/runWormCellTracking.sh '" 
                + matlabDirName +"' " + nGroups + " 0 1")
        commandList.insert(len(commandList)-1, qsubCommand1)
        print(qsubCommand1)
        f.write(qsubCommand1)
        f.write('\n')

        qsubCommand2 = ("qsub" + " -m as -l " + queueUsed + " -l job_mem_gb=" + "20" + " -l h_data=" 
                + "20" + "g -cwd" + " -o \"" + folderName + "\" -e \"" 
                + folderName + "\" -t 1 "
                + " -N \"" + jobName2 + "\" -hold_jid \"" + jobName + "\" "
                + " ../scripts/shae-pythonSubmissionScripts/runWormTrackCompiler.sh '" 
                + matlabDirName +"' '" + matlabDirName2 + "'")

        qsubCommand3 = ("qsub" + " -l " + queueUsed + " -l gb=" + memReq + " -l h_data=" 
                + memReq + "g -cwd" + " -o \"" + folderName + "\" -e \" " 
                + folderName + "\" -t 1-" + totalNeuronRuns + ":" + str(1) 
                + " -N \"" + jobName3 + "\" -hold_jid \"" + jobName2 + "\" "
                + " ../scripts/shae-pythonSubmissionScripts/runWormBotChecker.sh '" 
                + matlabDirName2 +"' " + groupsPerCell + " 0")

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
    
# bind enter key and button
master.b = Button(master, text="Enter", width=10, command=lambda:callback1(master=master))
master.b.grid(row=16,columnspan=2, sticky=W+E)
master.bind("<Return>", callback1b)

master.e1.focus_set()

mainloop()
