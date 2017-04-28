#!/usr/bin/python
import pickle
import slurmInput as slurm
import socket
import guiHelper as gu
import os


def make_gui():
    # need to make the master Tk window at the very beginning
    
    # load pickle file and get default values
    prevUser=gu.pickle_load()
    defaultName = prevUser['username']
    defaultDate = prevUser['date']
    defaultFolder = prevUser['folderName']
             
    master = gu.submitTK(rows=11,cols=2)
    
    #add each row with text inputs
    master.addGuiField("User Name",'user_name',defaultName)
    master.addGuiField("Parent Path",'parent_path','/tigress/LEIFER/PanNeuronal')
    master.addGuiField("Date of data",'date',defaultDate)
    master.addGuiField("DataFolderName",'folder_name',defaultFolder)
    master.addGuiCheck("Email",'email_flag',1)
    #make Enter button, tie it to the callback1
    master.addGuiButton("Enter",b_command=lambda:callback1(master=master))

    if  socket.gethostname()=='tigressdata.princeton.edu':
        master.addGuiButton("Select Folder",b_command=lambda:gu.selectFolder(master=master))
        
    return master
  
def submitScript(master=None):
    # get inputs from master
    username = master.e['user_name'].get()
    beginOfPath=master.e['parent_path'].get()
    date=master.e['date'].get()
    folderName=master.e['folder_name'].get()
    emailFlag       = master.e['email_flag'].var.get()

    # which folder to process, must add paths linux style
    fullPath = beginOfPath + "/" + date
    fullPath = fullPath + "/" + folderName
        
    print("Username: " + username)
    print("full path is: " + fullPath)

        
    # connect to della
    if 'password' in master.e.keys():
        password =master.e['password'].get
        client=gu.dellaConnect(username,password)
    else:
        client=gu.dellaConnect(username)
        
    #make a folder for the output files to be saved in
    slurm.make_ouputfolder(client,fullPath)

        
    #Build command list to be submitted to della
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    # add commands to the command list
    commandList=slurm.get_git_hash(commandList,client)
    commandList=slurm.path_setup(commandList)
    
    #submit flashfinder job
    commandList=slurm.flash_input(commandList,fullPath,emailFlag)
    
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
        Output files will be saved in 
        '''
        + fullPath 
        + ''' 
        and the corresponding LowMag folder
        The output files are:
                  hiResData.mat
                  cam0flashTrack.mat
                  cam1flashTrack.mat
                  ''')
                  
    #save defaults using pickle dump
    pickle_path = (os.environ['HOME'] + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    prevUser=gu.pickle_load()
    prevUser['username']=username
    prevUser['date'] = date
    prevUser['folderName']=folderName
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
                  
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
        This is the submission script for running time alignment and detcting flashes in videos. 
        The lowMag folders must be inside the corresponding BrainScanner folder on della. 
        The videos must each contain at least 2 flashes.
        
        
        For a quick test, run this code as follows:
        User Name: <your username>
        Parent Path:/tigress/LEIFER/PanNeuronal
        Date of Data: testing_sets
        Data Folder Name: Brain_working_dataset
        
        ''')
    master=make_gui()
    master.e['user_name'].focus_set()
    master.run()
