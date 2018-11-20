# submitCustomMatlab
# Run to submit any matlab function to run on della. A single job request will be made, so parfors can be used for multiple threads. The Code Input string will be exactly what is run in matlab. 



#!/usr/bin/python
import slurmInput as slurm
import socket
import guiHelper as gu
import os
import getpass


def make_gui():
    # need to make the master Tk window at the very beginning
    
    # load pickle file and get default values
    prevUser=gu.pickle_load()
    defaultName = prevUser['username']
    defaultTime = prevUser['time']
    defaultMemory = prevUser['mem']
    defaultDate = prevUser['date']
    defaultFolder = prevUser['folderName']
    defaultCommand=prevUser['matlab_command']
    defaultPath=prevUser['code_path']
    master = gu.submitTK(rows=11,cols=2)
    
    #add each row with text inputs
    master.addGuiField("User Name",'user_name',defaultName)
    
    #added new memory and time allotment
    master.addGuiField("Minutes Allotted",'time',defaultTime)
    master.addGuiField("mB Requested",'mem',defaultMemory)
    
    #select a .path file with all the folders to add to the matlab path. Each folder must end in /, and there must be an additional character return after the last entry. 
    master.addGuiField("Code Path",'code_path',defaultPath)
    master.addGuiField("Matlab Command",'matlab_command',defaultCommand)
    master.addGuiCheck("Email",'email_flag',1)
    #make Enter button, tie it to the callback1
    master.addGuiButton("Enter",b_command=lambda:callback1(master=master))

    if  socket.gethostname()=='tigressdata.princeton.edu':
        master.addGuiButton("Select Code Path File",b_command=lambda:gu.selectPath(master=master))
        
    return master
  
  
def submitScript(master=None):
    # get inputs from master
    username = master.e['user_name'].get()
    time=master.e['time'].get()
    mem=master.e['mem'].get()
    matlab_command=master.e['matlab_command'].get()
    emailFlag       = master.e['email_flag'].var.get()
    code_path= master.e['code_path'].get()
    
    #The user that is currently logged in, not the user that will be used to submit the job
    real_user= getpass.getuser()
    fullPath='/tigress/LEIFER/'+real_user
    
    print("Username: " + username)

        
    # connect to della
    if 'password' in master.e.keys():
        password =master.e['password'].get
        client=gu.dellaConnect(username,password)
    else:
        client=gu.dellaConnect(username)
        
    slurm.make_ouputfolder(client,fullPath)
    
    #Build command list to be submitted to della
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    #need to set path file for this

    # add commands to the command list
    commandList=slurm.get_git_hash(commandList,client)
    commandList=slurm.path_setup(commandList)
    
    pathCommand='export PATH_FILE='+ code_path
    commandList.insert(len(commandList)-1, pathCommand)  
      
      
    #submit custom input job
    commandList=slurm.custom_input(commandList,matlab_command,fullPath,emailFlag,time=time,mem=mem)
    
    #save defaults using pickle dump
    master.pickleDump()
    
    #write commands to text file via paramiko
    slurm.write_input(commandList,client,fullPath)

#print stdout and stderr to command line
    commands = "\n".join(commandList)
    stdin, stdout, stderr = client.exec_command(commands)
    print('stdOutput:')
    returnedOutput = stdout.readlines()
    print(' '.join(returnedOutput))
    print('stdError:')
    print(stderr.readlines())
    print('Done submitting job.\n\n')
# show user ooutput location
    print('''
                  ''')
                  
                  
    # close window at the end
    client.close()
    master.destroy()
    
        
def callback1(event=None,master=None):
    
    #Check for password and continue to submit job
    print(master.e['user_name'] .get())
    username = master.e['user_name'].get()
    isNeedsPassword=gu.passwordReqCheck(username)
    
    if isNeedsPassword:
        # use the same window as before, just add an additional password field
        # password
        master.addGuiField('password','password',default='******',show='*')
        master.addGuiButton("Enter",lambda:submitScript(master=master))
    else:
        print("No password needed")
        submitScript(master)
        

if __name__ == '__main__':
# bind enter key and button
    print('''
        The script submits any matlab function to be run on Della.
        The Matlab Command string will be excatly typed into matlab
        and run. The required code dependencies must be listed in a
        .path file. That file must have a list of all directories that
        have relavent code for your function. Each directory must end 
        in a "/" and the the file must end in a character return. 
        
        You can use parfor loops in matlab to multithread processes. 
        You may have to check exactly how to request resources.
        
        ''')
    master=make_gui()
    master.e['user_name'].focus_set()
    master.run()
