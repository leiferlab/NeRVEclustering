#!/usr/bin/python
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


if os.name == 'posix':
    os.system('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' ''')
        
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
    prevUser['username']=username
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    

    print("full path is: " + fullPath)
        
    
    # connect and submit job via sbatch
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
    commandList.insert(len(commandList)-1, "export PATH=$PATH:/usr/local/emboss/bin:/sge/bin/lx-amd64:/usr/kerberos/bin")

   # matlabDirName = "\\ ".join(matlabDirName);
    print('Writing inputs line to text file')
    userEmail=username+"@princeton.edu"
    
    commandList=slurm.flash_input(commandList,fullPath)
        
    slurm.make_ouputfolder(client,fullPath)
    #write commands to text file via paramiko
    slurm.write_input(commandList,client,fullPath)

    commands = "\n".join(commandList)
   # stdin, stdout, stderr = client.exec_command(commands)
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
    password = master.e7.get()
    
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
        # try to see if a password is needed
        client = paramiko.SSHClient()
        client.load_system_host_keys()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        username = master.e['user_name'].get()
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
        
        master.e['password'] = Entry(master,show="*")
        master.e['password'].insert(0, "Password")
        master.e['password'].grid(row=7, column=1, sticky=W+E)
        
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

master.e['user_name'].focus_set()

mainloop()
