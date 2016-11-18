#!/usr/bin/python
import pickle
import os
import sys
import subprocess
import tarfile
import time
import glob
from Tkinter import *
import paramiko
import re
import tkFileDialog
import math

# need to make the master Tk window at the very beginning
master = Tk()

# process input arguments
fileInput = sys.argv;
try:
    fileInput = fileInput[1]
except IndexError:
    fileInput = tkFileDialog.askdirectory(parent=master)

# which folder to process
fullPath,fileName = os.path.split((fileInput));
beginOfPath,folderName = os.path.split((fullPath));

print(("Base path: " + beginOfPath));
print(("Folder name: " + folderName));
print(("File name: " + fileName));


# store where we start from
startingDir = os.getcwd();
print((time.asctime()+'\n'));

# get ready for pickled variables 
pickle_path = (os.environ['HOME'] + "/platypusTemp/");
pickle_file = pickle_path + "pickles2.p";

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
	defaultName = "USER";


# use the example of the calculator from http://zetcode.com/gui/tkinter/layout/ to help layout
master.title("Options")

master.columnconfigure(0, pad=3);
master.columnconfigure(1, pad=3);
master.columnconfigure(2, pad=3);
master.columnconfigure(3, pad=3);
master.columnconfigure(4, pad=3);
master.columnconfigure(5, pad=3);
master.columnconfigure(6, pad=3);

master.rowconfigure(0, pad=3);
master.rowconfigure(1, pad=3);

# user name
master.L1 = Label(master, text="User Name")
master.L1.grid(row=0, column=0, sticky=W+E)

master.e1 = Entry(master)
master.e1.insert(0, defaultName)
master.e1.grid(row=0, column=1, sticky=W+E)

# memory per cell
master.L2 = Label(master, text="Memory per cell (gb)")
master.L2.grid(row=1, column=0, sticky=W+E)

master.e2 = Entry(master)
master.e2.insert(0,'2')
master.e2.grid(row=1, column=1, sticky=W+E)

# grouping size
master.L3 = Label(master, text="Group Size")
master.L3.grid(row=2, column=0, sticky=W+E)

master.e3 = Entry(master)
master.e3.insert(0,'5')
master.e3.grid(row=2, column=1, sticky=W+E)

# use 1day queue
master.L4 = Label(master, text="Use 1 Day Queue")
master.L4.grid(row=3, column=0, sticky=W+E)

var = IntVar();
master.CB4 = Checkbutton(master, text=None, variable=var);
master.CB4.var = var;
master.CB4.grid(row=3, column=1, sticky=W+E);

# upload data
master.L5 = Label(master, text="update mat files")
master.L5.grid(row=4, column=0, sticky=W+E)

var2 = IntVar();
master.CB5 = Checkbutton(master, text=None, variable=var2);
master.CB5.var2 = var2;
master.CB5.grid(row=4, column=1, sticky=W+E);

# Number of Frames
master.L6 = Label(master, text="Number of Points")
master.L6.grid(row=5, column=0, sticky=W+E)

master.e6 = Entry(master)
master.e6.insert(0,'200')
master.e6.grid(row=5, column=1, sticky=W+E)

if os.name == 'posix':
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        
def submitScript(master=None):
    # username
    username = master.e1.get();
    print("Username: " + username);
    pickle.dump( { "username": username}, open(pickle_file, "wb" ) )
    
    # memory requested
    memReq = master.e2.get();
    print("Memory requested: " + memReq + " gb per cell.");
    
    # grouping size
    groupSize = master.e3.get();
    print("Group Size is" + groupSize);

    numberOfFrames=master.e6.get();
    print("number of Frames"+ numberOfFrames)
    # update code
    if master.CB4.var.get()==1:
        queueUsed = "1day";
        print("Request: 1 day queue")
    else:
        queueUsed = "1hr";
        print("Request: 1 hr queue")
        
    if master.CB5.var2.get()==1:
         uploadCells = 1;
         print("Request: Uploading data")
    else:
         uploadCells=0;
         print("Request: Not Uploading data")


    isUpdateCode = 0;
    
    if isUpdateCode == 1:
        # see stackoverflow http://stackoverflow.com/questions/13851846/paramiko-sftpclient-setting-missing-host-key-policy
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
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
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        try:
            client.connect('cetus.princeton.edu', 22, username)
        except paramiko.AuthenticationException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
        except SSHException:
            password = master.e7.get()
            client.connect('cetus.princeton.edu', 22, username, password)
      
        commandList = ["pwd","pwd"]; # pwd at both ends, give the list something to add to the middle of
        commandList.insert(len(commandList)-1, "export GRID_HOME=$HOME");
        commandList.insert(len(commandList)-1, "cd scripts/shae-pythonSubmissionScripts");
        commandList.insert(len(commandList)-1, "echo environment variables");
        commandList.insert(len(commandList)-1, "echo $GRID_HOME");
        commandList.insert(len(commandList)-1, "python updateCodeFromGit.py");
        commands = "\n".join(commandList);
        #print(commands)
        stdin, stdout, stderr = client.exec_command(commands)
        print('Updating code from Git')
        returnedOutput = stdout.readlines();
        print(' '.join(returnedOutput))
        print('stdError: updating from Git')
        print(stderr.readlines())
        client.close()
        print('Done updating code.\n\n')
        
        
    if uploadCells == 1: 
        os.chdir(beginOfPath);
        
        # copy images over via sftp
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        try:
         client.connect('cetus.princeton.edu', 22, username)
        except paramiko.AuthenticationException:
         password = master.e7.get()
         client.connect('cetus.princeton.edu', 22, username, password)
        except SSHException:
         password = master.e7.get()
         client.connect('cetus.princeton.edu', 22, username, password)
         
        client2 = client.open_sftp()
        client2.chdir('.'); # change the current director to itself, not sure why this is necessary
        # move the file, hopefully to the right location
        try:
            client2.stat('data/' + folderName)
        except:
            client2.mkdir('data/' + folderName)

        print(os.path.join(folderName,fileName));
        client2.put(os.path.join(fullPath,fileName),("data/"+os.path.join(folderName,fileName)))
                            
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
    matlabDirName = os.path.join(folderName,fileName);
    jobName = "_".join(matlabDirName);
    # matlabDirName = "\\ ".join(matlabDirName);
    print(matlabDirName)
    
    # connect and submit job via qsub
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.WarningPolicy)
    try:
        client.connect('cetus.princeton.edu', 22, username)
    except paramiko.AuthenticationException:
        password = master.e7.get()
        client.connect('cetus.princeton.edu', 22, username, password)
    except SSHException:
        password = master.e7.get()
        client.connect('cetus.princeton.edu', 22, username, password)
    
    commandList = ["pwd","pwd"]; # pwd at both ends, give the list something to add to the middle of
    # set up the environment so that it matches an ssh login instead of the reduced paramiko one, hopefully this will help.
    #commandList.insert(len(commandList)-1, "export MANPATH=/usr/share/man:/usr/local/man:/sge/man");
    #commandList.insert(len(commandList)-1, "export HOSTNAME=cetus.gen.pvt");
    #commandList.insert(len(commandList)-1, "export FAFNER_HOME=/Genomics/grid/users/$USER");
    #commandList.insert(len(commandList)-1, "export TERM=xterm-256color");
    #commandList.insert(len(commandList)-1, "export HISTSIZE=1000");
    #commandList.insert(len(commandList)-1, "export WUBLASTDB=/grid/data/blastdb");
    #commandList.insert(len(commandList)-1, "export SGE_ARCH=lx-amd64");
    #commandList.insert(len(commandList)-1, "export SGE_CELL=default");
    #commandList.insert(len(commandList)-1, "export BLASTFILTER=/usr/local/wu/filter");
    #commandList.insert(len(commandList)-1, "export SSH_TTY=/dev/pts/7");
    #commandList.insert(len(commandList)-1, "export LD_LIBRARY_PATH=/opt/gurobi55/linux64/lib");
    #commandList.insert(len(commandList)-1, "export LS_COLORS=");
    #commandList.insert(len(commandList)-1, "export _XX_SGE_BASHRC=1");
    #commandList.insert(len(commandList)-1, "export WUBLASTMAT=/usr/local/wu/matrix");
    #commandList.insert(len(commandList)-1, "export OS_ARCH=linux-amd64");
    #commandList.insert(len(commandList)-1, "export INPUTRC=/etc/inputrc");
    #commandList.insert(len(commandList)-1, "export LANG=en_US.UTF-8");
    commandList.insert(len(commandList)-1, "export SGE_ROOT=/sge");
    #commandList.insert(len(commandList)-1, "export GUROBI_HOME=/opt/gurobi55/linux64");
    #commandList.insert(len(commandList)-1, "export BLASTDB=/grid/data/blastdb");
    #commandList.insert(len(commandList)-1, "export WUBLASTFILTER=/usr/local/wu/filter");
    #commandList.insert(len(commandList)-1, "export EMBOSS_ACDROOT=/usr/local/emboss/share/EMBOSS/acd");
    #commandList.insert(len(commandList)-1, "export CVS_RSH=ssh");
    #commandList.insert(len(commandList)-1, "export LESSOPEN=|/usr/bin/lesspipe.sh %s");
    #commandList.insert(len(commandList)-1, "export BLASTMAT=/usr/local/wu/matrix");
    #commandList.insert(len(commandList)-1, "export G_BROKEN_FILENAMES=1");
    #commandList.insert(len(commandList)-1, "export MAIL=/var/spool/mail/$USER");
    commandList.insert(len(commandList)-1, "export PATH=$PATH:/usr/local/emboss/bin:/sge/bin/lx-amd64:/usr/kerberos/bin");
    commandList.insert(len(commandList)-1, "export GRID_HOME=$HOME");
    commandList.insert(len(commandList)-1, "cd data");
    os.chdir(beginOfPath);
    
    #deal with folders with spaces in names

    matlabDirName = os.path.join(folderName,fileName);
    jobName = folderName;
   # matlabDirName = "\\ ".join(matlabDirName);
    print(folderName);
    qsubCommand = ("qsub" + " -m as -l " + queueUsed + " -l gb=" + memReq + " -l h_data=" 
		+ memReq + "g -cwd" + " -o \"" + folderName + "\" -e \"" 
		+ folderName + "\" -t 1-" + numberOfFrames + ":" + groupSize 
		+ " -N \"" + jobName + "\" ../scripts/shae-pythonSubmissionScripts/runWormBotChecker.sh '" 
		+ matlabDirName +"'");	
    print(qsubCommand);
    commandList.insert(len(commandList)-1, qsubCommand);
    commands = "\n".join(commandList);
    stdin, stdout, stderr = client.exec_command(commands)
    print('stdOutput: submitting job')
    returnedOutput = stdout.readlines();
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
            master = event.widget.master;
        
        username = master.e1.get();
        print(username);
        # password, just the number of characters
        password = master.e7.get();
        
        passList = ["*"] * len(password);
        print("Entered password: " +   "".join(passList));
        
        submitScript(master)
    
    username = master.e1.get();
    if username == "noPass":
        isNeedsPassword = True;
    else:
        # try to see if a password is needed
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.WarningPolicy)
        username = master.e1.get()
        try:
            client.connect('cetus.princeton.edu', 22, username)
            isNeedsPassword = False
        except paramiko.AuthenticationException:
            isNeedsPassword = True;
            
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
    master = event.widget.master;
    #master.withdraw()
    callback1(master=master)
    
# bind enter key and button
master.b = Button(master, text="Enter", width=10, command=lambda:callback1(master=master))
master.b.grid(row=6,columnspan=2, sticky=W+E)
master.bind("<Return>", callback1b)

master.e1.focus_set()

mainloop()
