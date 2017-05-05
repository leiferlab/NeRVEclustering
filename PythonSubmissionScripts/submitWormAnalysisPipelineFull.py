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
    master.addGuiField("Number of frames",'nframes',defaultFrameNumber)
    master.addGuiField("Number of Reference",'n_ref',defaultRefNumber)
    master.addGuiField("Number of Neurons",'n_neurons',defaultNeuronNumber)
    master.addGuiField("Number of checks",'n_checks',defaultCheckNumber)
    # add check box inputs
    master.addGuiCheck("Run Straightening",'straight_flag')
    master.addGuiCheck("Run NeRVE",'track_flag')
    master.addGuiCheck("Run Error Correction",'check_flag')
    master.addGuiCheck("Run Signal Extraction",'crop_flag')
    master.addGuiCheck("Send Email",'email_flag',1)
    #make Enter button, tie it to the callback1
    master.addGuiButton("Enter",b_command=lambda:callback1(master=master))

    #if we're logged onto tigressdata, then add a select folder button to navigate files
    if  socket.gethostname()=='tigressdata.princeton.edu':
        print("Navigate inside a Brainscanner folder and press Select.")
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
        output_string+=['Straightening Result: CLstraight folder']
        
    if trackFlag:
        commandList=slurm.track_input(commandList,fullPath,n_volumes,n_references)
        output_string+=['NeRVE Result: trackMatrix folder and pointsStats2.mat']

    if checkFlag:
        commandList=slurm.check_input(commandList,fullPath,n_volumes,n_check,n_neurons)
        output_string+=['Error Checking Result: BotCheckFolder and pointsStatsNew.mat']

    if cropFlag:
        commandList=slurm.crop_input(commandList,fullPath,emailFlag)
        output_string+=['Signal Extractoin Result: heatData.mat and botFiducials Folder']
        
        
    #save inputs using pickle dump, they will appear as defaults at the next call. 
    master.pickleDump()
    
    #write command list to text file via paramiko
    slurm.write_input(commandList,client,fullPath)

    #submit the command list to della
    commands = "\n".join(commandList)
    stdin, stdout, stderr = client.exec_command(commands)
    returnedOutput = stdout.readlines()
    returnedErr = stderr.readlines()
    
    returnedOutput=['#### OUTPUT \n']+returnedOutput
    returnedErr=['#### ERROR \n']+returnedErr
    
    slurm.write_input(returnedOutput,client,fullPath)
    slurm.write_input(returnedErr,client,fullPath)


    print('stdOutput:')
    print(' '.join(returnedOutput))
    print('stdError:')
    print(' '.join(returnedErr))
    
    print('Done submitting job.\n\n')
    
    
    
    print('''
        The following output files will be produced and saved in
        '''
        + fullPath
        + ':\n'
        + '\n'.join(output_string))
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
        

    
    
# actually runs the code
if __name__ == '__main__':
# bind enter key and button
    print('''
        This is the submission script for running analysis on whole brain imaging.
        This code runs straightening, tracking, and cross validation of points. It also 
        extracts signals and produces heatmaps. This can be run after the wormAnalysisPreview
        gui gives all green lights. Check the readme for instructions for della access. 
        
        Requirements:
                centerline.mat: created by submitWormAnalysisCenterline.py
                    must be located in the behavior folder
                alignments.mat: from alignment_gui, located in Brainscanner folder
                all timing results from submitWormFlashFinder.py
                
                
        For a quick test, run this code as follows:
        User Name: <your username>
        Parent Path:/tigress/LEIFER/PanNeuronal
        Date of Data: testing_sets
        Data Folder Name: Brain_working_dataset
        Number of Frames : All  (or specify some subset of frames for testing.. e.g. 1000)
        Number of References: 10 (This determines how good our tracking is. 300 was used by default. For stationary animals, use ~25)
        Number of Neurons: 150 (Sets the upper limit for the number of neurons in the worm)
        Number of Checks: 100 (Number of volumes to use for the error checking step... ~100-300 is normally fine)
        <click all check boxes>  (you would only uncheck the Run Straighting, Run Track, Run Check or Run Crop for debugging purposes)

        
        ''')
    master=make_gui()
    master.e['user_name'].focus_set()
    master.run()
