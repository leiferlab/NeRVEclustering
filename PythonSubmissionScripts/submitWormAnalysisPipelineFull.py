#!/usr/bin/python



# code for submitting jobs onto Princeton Slurm cluster,  using Della now.
# Running this opens a GUI for selecting input files and parameters for tracking
# Neurons in freely moving c. elegans videos. The neuron tracking algorithm
# is outlined in Nguyen et all (2016). This is the main part of the algorithm, with
# straightening and neuron registration vector encoding. It is currently best to run
# this program from terminal on tigressdata. The matlab code that is called saves
# analysis output files into the same folder specified by the GUI. The final output
# is the heatData.mat file, which has behavior and fluorescence data for the worm. 

# This code runs after centerlines are found (submitWormAnalysisCenterlines.py) and 
# timings are found (submitWromFlashFinder.py). 

# Jeffrey Nguyen

import pickle
import slurmInput as slurm
import socket
import guiHelper as gu
import os


# need to make the master Tk window at the very beginning
def make_gui():
    
    # load pickle defaults
    prevUser=gu.pickle_load()
    defaultName = prevUser['username']
    defaultDate = prevUser['date']
    defaultFolder = prevUser['folderName']
    defaultFrameNumber = prevUser['frameNumber']
    defaultRefNumber = prevUser['refNumber']
    defaultNeuronNumber = prevUser['neuronNumber']
    defaultCheckNumber = prevUser['checkNumber']

    #build the initial GUI
    master = gu.submitTK(rows=19,cols=2)
    
    #add each row with text inputs
    master.addGuiField("User Name",'user_name',defaultName)
    master.addGuiField("Parent Path",'parent_path','/tigress/LEIFER/PanNeuronal')
    master.addGuiField("Date of data",'date',defaultDate)
    master.addGuiField("DataFolderName",'folder_name',defaultFolder)
    master.addGuiField("number of frames",'nframes',defaultFrameNumber)
    master.addGuiField("Number of Reference",'n_ref',defaultRefNumber)
    master.addGuiField("Number of Neurons",'n_neurons',defaultNeuronNumber)
    master.addGuiField("Number of checks",'n_checks',defaultCheckNumber)
    # add check box inputs
    master.addGuiCheck("Run Straightening",'straight_flag')
    master.addGuiCheck("Run Track",'track_flag')
    master.addGuiCheck("Run Check",'check_flag')
    master.addGuiCheck("Crop",'crop_flag')
    master.addGuiCheck("Email",'email_flag',1)
    #make Enter button, tie it to the callback1
    master.addGuiButton("Enter",b_command=lambda:callback1(master=master))

    if  socket.gethostname()=='tigressdata.princeton.edu':
        master.addGuiButton("Select Folder",b_command=lambda:gu.selectFolder(master=master))
    return master


# runs on enter button, gets gui inputs, connects to della, and submits jobs
def submitScript(master=None):
    
    # get inputs from master
    username = master.e['user_name'].get()
    beginOfPath=master.e['parent_path'].get()
    date=master.e['date'].get()
    folderName=master.e['folder_name'].get()
    n_references=master.e['n_ref'].get()
    n_neurons=master.e['n_neurons'].get()
    n_check=master.e['n_checks'].get()

    # construct path to data folder
    fullPath = beginOfPath + "/" + date + "/" + folderName
    
    #Which parts of code to run
    straightFlag    = master.e['straight_flag'].var.get()
    trackFlag       = master.e['track_flag'].var.get()
    checkFlag       = master.e['check_flag'].var.get()
    cropFlag        = master.e['crop_flag'].var.get()
    emailFlag       = master.e['email_flag'].var.get()
    
    # connect to della
    if 'password' in master.e.keys():
        password =master.e['password'].get
        client=gu.dellaConnect(username,password)
    else:
        client=gu.dellaConnect(username)
        
        
    #if n_volumes is a string (normall "all"), read file from submission_parameters.txt in the Brain scannner folder.
    # otherwise, just get the value
    n_volumes=master.e['nframes'].get()
    if not n_volumes.isdigit():
        n_volumes=gu.get_nframes(client,fullPath)
    n_volumes=int(n_volumes)

    print("Username: " + username)    
    print("Number of Neruons is " + n_neurons)
    print("Total runs is " + str(n_volumes))

    #Construct command list
    commandList = ["pwd","pwd"] # pwd at both ends, give the list something to add to the middle of
    # set up the environment so that it matches an ssh login instead of the reduced paramiko one, hopefully this will help.
    slurm.make_ouputfolder(client,fullPath)
        
    #submit path setup bash commands and add display of git hash 
    commandList=slurm.get_git_hash(commandList,client)
    commandList=slurm.path_setup(commandList)
    
    # also make output string to tell user what files to expect to see
    output_string=[]
    #add sbatch command to commandList for each flag.
    if straightFlag:
        commandList=slurm.straighten_input(commandList,fullPath,n_volumes)
        output_string+='CLstraight folder'
        
    if trackFlag:
        commandList=slurm.track_input(commandList,fullPath,n_volumes,n_references)
        output_string+='trackMatrix folder'
        output_string+='pointsStats2.mat'

    if checkFlag:
        commandList=slurm.check_input(commandList,fullPath,n_volumes,n_check,n_neurons)
        output_string+='trackMatrix folder'
        output_string+='pointsStatsNew.mat'

    if cropFlag:
        commandList=slurm.crop_input(commandList,fullPath,emailFlag)
        output_string+='heatData.mat'
        output_string+='botFiducials Folder'

    #submit the command list to della
    commands = "\n".join(commandList)
    stdin, stdout, stderr = client.exec_command(commands)
    #write commands to text file via paramiko
    slurm.write_input(commandList,client,fullPath)
    
    print('stdOutput:')
    returnedOutput = stdout.readlines()
    print(' '.join(returnedOutput))
    print('stdError:')
    print(stderr.readlines())
    print('Done submitting job.\n\n')
    
    print('''
        Output files will be saved in 
        '''
        + fullPath
        + '\n'
        + '\n'.join(output_string))
    # close window at the end
    
    #save inputs using pickle dump, they will appear as defaults at the next call. 
    pickle_path = (os.path.expanduser('~') + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    prevUser=gu.pickle_load()
    prevUser['username']=username
    prevUser['frameNumber']=n_volumes
    prevUser['date'] = date
    prevUser['folderName']=folderName
    prevUser['refNumber']=n_references
    prevUser['neuronNumber']=n_neurons
    prevUser['checkNumber']=n_check
    pickle.dump(prevUser, open(pickle_file, "wb" ) )
    
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
        This is the submission script for running analysis on whole brain imaging.
        This code runs straightening, tracking, and cross validation of points. It also 
        extracts signals and produces heatmaps. This can be run after the wormAnalysisPreview
        gui gives all green lights. Check the readme for instructions for della access. 
        
        
        For a quick test, run this code as follows:
        User Name: <your username>
        Parent Path:/tigress/LEIFER/PanNeuronal
        Date of Data: testing_sets
        Data Folder Name: BrainScanner20161031_111303
        Number of Frames : All
        Number of References: 10
        Number of Neurons: 150
        Number of Checks: 100
        <click all check boxes>
        
        ''')
    master=make_gui()
    master.e['user_name'].focus_set()
    master.run()
