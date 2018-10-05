#!/usr/bin/python
# Supplementary code to go with submitWormAnalysis code developed in the Leifer Lab. For the most part, the code is generating sbatch jobs with flags and inputs. For instructions for the flags, see
# https://slurm.schedmd.com/sbatch.html
# but the general strcutre will go as
# sbatch 
#    --mem=<memory request> 
#    --time = < time request> 
#    -D < directory to start>
#    -J < name of job>
#    -d < dependency>
#    --output =< output file path>
#    --error = < error file path>
#    < email string>
#    --array = 1-<end>:stepsize
#    < shell script that runs job>
#    < inputs for the shell script>
# 
#  functions here accompany submitWormAnalysis codes 
# Jeffrey Nguyen

import os
import numpy as np
import datetime
import pickle
import subprocess
import getpass
import socket
import time

CODE_PATH='/tigress/LEIFER/communalCode/3dbrain' #path to the code repo
MIN_TIME_STR = "--time=6:02:00"  #minimum time string for use on short queue
PS_NAME1 =  'PointsStats.mat'
PS_NAME2 =  'PointsStats2.mat'
NOW=datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") # datetime string

# apparently an update of slurm and/or della requires now an absolute path after 
# "-D" in slurm submission commands. Replacing os.path.basename(fullPath) with a
# call to this function, so that it can be reverted easily.
def get_folder_name(fullPath):
    # it was
    #folderName=os.path.basename(fullPath)
    folderName=fullPath
    return folderName

# construct email string using the user currently logged on. This is fine if run from tigressdata, but may have problems when run from home computers where the user is not a princeton netID. 
def get_email_script(mail_type='end,fail'):
    user= getpass.getuser()
    user_email=user+'@princeton.edu'
    mail_script= ' --mail-type=' + mail_type +' --mail-user=' + user_email
    return mail_script

# make a path name for the output file where the .out and .err files will be saved for all the jobs. One is created each day and all files from that day are saved in the folder. 
def make_output_path(fullPath):
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    return outputFilePath
    

# submit command over ssh to get the current git has of CODE_PATH, add it to the commandList so that it will appear in the input.txt file
def get_git_hash(commandList,client):
    git_command='cd '+ CODE_PATH +'\n git rev-parse HEAD \n git branch -vv \n cd $HOME'
    stdin,stdout,stderr = client.exec_command(git_command)
    current_git_hash=stdout.readlines()
    commandList.insert(len(commandList)-1, '#### current git: '+current_git_hash[0])
    commandList.insert(len(commandList)-1, '#### current git: '+current_git_hash[1])
    
    return commandList
    
# set up path so .sh files know where to find the matlab codes
def path_setup(commandList):
    # add a header
    commandList.insert(len(commandList)-1, '####PATH SETUP####'+NOW)
    code_home,_=os.path.split(CODE_PATH)
    # export CODE_HOME, needed for the shell scripts to add the correct path for matlab
    commandList.insert(len(commandList)-1, "export CODE_HOME="+code_home)
    commandList.insert(len(commandList)-1, "umask 002") #make permissions open to group and user for all future files
    return commandList
    
    
# input code for initializing the centerline worksapce
def centerline_start_input(commandList,fullPath,email_flag=False):
    commandList.insert(len(commandList)-1, '####CENTERLINES START####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath= make_output_path(fullPath)
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    
    input0 = "clusterCL_start('"+ fullPath + "')"
    #make submission command
    qsubCommand0 = ("sbatch --mem=18000 " 
        + "--time=2:02:00"
        + " -N 1 -n 1 -c 10" #for use of 10 cores for parfor loops in the matlab code
        + " -D " + folderName
        + " -J "+ folderName
        + " " + get_email_script('fail')
        + " --output=\"" + outputFilePath + "/CLstart-%J.out"+"\""
        + " --error=\"" + outputFilePath + "/CLstart-%J.err" + "\""
        + " " + code_runinput
        + " \"" + input0 +"\" ")
        
    commandList.insert(len(commandList)-1, qsubCommand0)
    return commandList
        
# input code for centerline submission
def centerline_input(commandList,fullPath,email_flag = False,chip_flag = 0):
    #make header for input file
    commandList.insert(len(commandList)-1, '####CENTERLINES####'+NOW)
    
    folderName=get_folder_name(fullPath)
    outputFilePath= make_output_path(fullPath)
    
    #path to shell script job files
    code_centerline = CODE_PATH + '/PythonSubmissionScripts/runWormCenterlineFitting.sh'
    code_centerline_compile = CODE_PATH + '/PythonSubmissionScripts/runWormCenterlineCompile.sh'
    
    if email_flag:
        email_script=get_email_script('end,fail')
    else:
        email_script=get_email_script('fail')
        
    
    qsubCommand1 = ("sbatch --mem=18000 " 
        + "--time=2:02:00" 
#        + "--time=7320"
        + " -D " + folderName
        + " -J "+ folderName
        + " -d singleton" #dependency (only allow one job with this name at a time.. in this case wait for the previous job to complete)
        + " --output=\"" + outputFilePath + "/CLjob-%J.out"+"\"" #where output messages are stored
        + " --error=\"" + outputFilePath + "/CLjob-%J.err" + "\"" #where error messages are stored
        + " --array=1-32:1"  #there are 32 threads for centerline fitting
        + " " + code_centerline #this is the shell command to invoke the centerline specific matlab script
        + " '"  + fullPath +"' "+ str(chip_flag) +" ") #this is an input which is handed into the shell command which gets into matlab
        
    commandList.insert(len(commandList)-1, qsubCommand1) #Add this to the stack of strings
   
    #Same as above but now for the "centerline compile" part 
    qsubCommand2 = ("sbatch --mem=18000 " 
#        + MIN_TIME_STR 
        + "--time=2:02:00"
        + " -D " + folderName
        + " -J "+ folderName
        + email_script
        + " -d singleton"
        + " --output=\"" + outputFilePath + "/CLCompile-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/CLCompile-%J.err" + "\""
        + " " + code_centerline_compile
        + " '"  + fullPath +"' ")
        
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    return commandList


def straighten_input(commandList,fullPath,totalRuns,email_flag = False):
    commandList.insert(len(commandList)-1, '####STRAIGHTENING####'+NOW)

    totalRuns = int(totalRuns)
    folderName=get_folder_name(fullPath)
    outputFilePath= make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'
    code_straighten = CODE_PATH + '/PythonSubmissionScripts/runWormStraighten.sh'
    code_pscompiler = CODE_PATH + '/PythonSubmissionScripts/runWormCompilePointStats.sh'
    fail_script=email_script=get_email_script('fail')

    if email_flag:
        email_script=get_email_script()
    else:
        email_script=get_email_script('fail')
    
    input0 = "clusterStraightenStart('"+ fullPath + "')"
    qsubCommand0 = ("sbatch --mem=16000 " 
        + MIN_TIME_STR 
        + " -D " + folderName
        + " -J "+ folderName
        + fail_script
        + " --output=\"" + outputFilePath + "/straight_s-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/straight_s-%J.err" + "\""
        + " " + code_runinput
        +" \"" + input0 +"\"")
    qsubCommand0 = "q0=$("+ qsubCommand0 + ")" # wrap the whole comaand in parens and set that equal to q0 at the bash level.
    commandList.insert(len(commandList)-1, qsubCommand0) #In the command list, have a command that runs the job and coppies the name into $q0 so that later we can submit jobs that depend on it and refer to it by its number stored in q0
    commandList.insert(len(commandList)-1, "echo $q0")
    
    dependencyString=" --dependency=afterok:${q0##* }" #wait for q0 (previous command) to be done successfully. The next job will use the dependency string.


#since we can only submit up to 1000 jobs, we have to submit nRuns number of seperate sbatch jobs
    nRuns=totalRuns//1000+1
    #arbitrary stepsize, these jobs are fairly fast, but for counting, limit to 300 jobs
    stepSize=totalRuns//300+1
    
    
    for i in range(0,int(nRuns)):
        #add offset as input to shell script so matlab can be called normally
        offset=str(i*1000)
        #last term for each sbatch, will be 1000 for all but the last sbatch command
        if i>=(nRuns-1):
            currentLimit=str(int(totalRuns)%1000)
        else:
            currentLimit="1000"

        qsubCommand1 = ("sbatch --mem=16000 " 
            + MIN_TIME_STR 
            + " -D " + folderName
            + " -J "+ folderName 
            + dependencyString
            + fail_script
            + " --output=\"" + outputFilePath + "/straight-%J.out" + "\""
            + " --error=\"" + outputFilePath + "/straight-%J.err" + "\""
            + " --array=1-" + currentLimit + ":" + str(stepSize) 
            + " " + code_straighten 
            + " '"  + fullPath +"' "  + str(stepSize) +" " + offset)
        commandList.insert(len(commandList)-1, qsubCommand1)
        
    qsubCommand2 = ("sbatch --mem=12000 " 
        + MIN_TIME_STR 
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + " " + email_script
        + " --output=\"" + outputFilePath + "/pscompile-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/pscompile-%J.err" + "\""
        + " " + code_pscompiler 
        +" '" + fullPath +"'")
    commandList.insert(len(commandList)-1, qsubCommand2)
    return commandList


#This one looks different because it needs to work around limits on thread names imposed by SLURM/ Della.
# Basically every thread needs a number. SLURM/Della forbits numbers greater than 1000 for thread names. 
# Because we want to name our threads in a way that it is obvious which volumes they handle, we often have to exceed the 1000 thread name limit. So in that case we run multiple submissions and get additional digits that way.  
def track_input(commandList,fullPath,totalRuns,nRef,email_flag = False):
    commandList.insert(len(commandList)-1, '####TRACKING####'+NOW)
    totalRuns=int(totalRuns)
    nRef=int(nRef)
    
    folderName=get_folder_name(fullPath)
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    
    if email_flag:
        email_script=get_email_script()
    else:
        email_script=get_email_script('fail')
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    code_track = CODE_PATH + '/PythonSubmissionScripts/runWormCellTracking.sh'
    code_trackcompiler = CODE_PATH+'/PythonSubmissionScripts/runWormTrackCompiler.sh'
    
    #currently arbitrary time request, ask for 2*nRef minutes (always above 180 min) 
    qString_track = "--time=" + str(np.max((nRef*2,180)))     

    matlabDirName = fullPath + "/" +  PS_NAME1
    matlabDirName2 = fullPath + "/" + PS_NAME2
    
    nRuns=totalRuns//1000+1
    stepSize=nRuns #Calculate how many volumes each parallel thread should handle serially. This is currently capped so that a maximium of 1500 threads are run in this section. 
    input1= "makePointStatsRef('"+ fullPath +"',"+ str(nRef) + ")"
    qsubCommand0 = ("sbatch --mem=2000 " 
        + MIN_TIME_STR 
        + " -D " + folderName
        + " -J "+ folderName
        + " -d singleton"
        + " --output=\"" + outputFilePath + "/track_s-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/track_s-%J.err" + "\""
        + " " + code_runinput
        +" \"" + input1 +"\"")
    qsubCommand0 = "q1=$("+ qsubCommand0 + ")"
    commandList.insert(len(commandList)-1, qsubCommand0)
    commandList.insert(len(commandList)-1, "echo $q1")
    commandList.insert(len(commandList)-1, '\r')
    #get jobID number to set up dependency for next job
    dependencyString=" --dependency=afterok:${q1##* }"
    

    for i in range(nRuns):
        #number of jobs capped at 1000, work around by submitting many batches with length 1000. 
        offset=str(i*1000)
        if i>=(nRuns-1):
            currentLimit=str(totalRuns%1000)
        else:
            currentLimit="1000"

        qsubCommand1 = ("sbatch --mem=4000 " 
            + qString_track 
            + dependencyString
            + " -J "+ folderName
            + " --output=\"" + outputFilePath + "/track-%J.out" + "\"" 
            + " --error=\"" + outputFilePath + "/track-%J.err" + "\""
            + " --array=1-" + currentLimit + ":" + str(stepSize)
            + " " + code_track + " '" 
            + fullPath +"' " + offset + " " + str(stepSize))
        
        commandList.insert(len(commandList)-1, qsubCommand1)
    commandList.insert(len(commandList)-1, '\r')
    
    qsubCommand2 = ("sbatch --mem=128000 " 
        + qString_track 
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + email_script
        + " --output=\"" + outputFilePath + "/trackCompiler-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/trackCompiler-%J.err" + "\""
        + " " + code_trackcompiler 
        +" '" + matlabDirName +"' '" + matlabDirName2 + "'")
    commandList.insert(len(commandList)-1, qsubCommand2)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


def check_input(commandList,fullPath,totalRuns,nCheck,nNeurons,email_flag = False):
    commandList.insert(len(commandList)-1, '####CHECKING####'+NOW)
    nCheck=int(nCheck)
    nNeurons=int(nNeurons)
    
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    if email_flag:
        email_script=get_email_script('end,fail')
    else:
        email_script=""
        
    code_checkcompiler = CODE_PATH + '/PythonSubmissionScripts/runWormBotCheckCompiler.sh'
    code_check = CODE_PATH+'/PythonSubmissionScripts/runWormBotChecker.sh'
    time_estimate=int(np.round(nCheck*totalRuns*.2/60))
    qString_check = "--time="+ str(np.max((time_estimate,180)))   

    matlabDirName2 = fullPath + "/" + PS_NAME2
    
    qsubCommand5 = ("sbatch --mem=8000 " 
        + qString_check 
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + " --output=\"" + outputFilePath + "/check-%J.out" + "\" "
        + " --error=\"" + outputFilePath + "/check-%J.err" + "\" "
        + " --array=1-"+ str(nNeurons) + ":" + "1"
        + " " + code_check 
        + " '"+ matlabDirName2 +"' "
        + str(nCheck))
    commandList.insert(len(commandList)-1, qsubCommand5)
        
    qsubCommand6 = ("sbatch --mem=16000 " 
        + "--time=10:02:00"
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + email_script
        + " --output=\"" + outputFilePath + "/checkcompile-%J.out" + " \""
        + " --error=\"" + outputFilePath + "/checkcompile-%J.err" + " \""
        + " " + code_checkcompiler 
        + "  '" + fullPath +"' ")
        
    commandList.insert(len(commandList)-1, qsubCommand6)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


def crop_input(commandList,fullPath, email_flag = False):
    commandList.insert(len(commandList)-1, '####CROPPING####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    
    input1= "fiducialCropper3('"+ fullPath +"')"

    if email_flag:
        email_script=get_email_script('end,fail')
    else:
        email_script=""
        
    qsubCommand7 = ("sbatch --mem=8000 " 
        + MIN_TIME_STR 
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + email_script
        + " --output=\"" + outputFilePath + "/crop-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/crop-%J.err" + "\""
        + " " + code_runinput 
        + " \"" + input1 +"\" ")
    
    commandList.insert(len(commandList)-1, qsubCommand7)
    commandList.insert(len(commandList)-1, '\r')
    return commandList



def flash_input(commandList,fullPath, email_flag = False):
    commandList.insert(len(commandList)-1, '####TIME SYNC####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'

    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""

    input1= "highResTimeTraceAnalysisTriangle4('"+ fullPath + "')"
    input2= "multipleAVIFlash('"+ fullPath +"')"
    qsubCommand1 = ("sbatch --mem=4000 " 
        + "--time=3:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/datFlash-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/datFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input1 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=2000 "  
        + "--time=2:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/avFlash-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/avFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input2 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    return commandList

# For new file format from LabView's new software
def flashNew_input(commandList,fullPath, email_flag = False):
    commandList.insert(len(commandList)-1, '####TIME SYNC####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'
    code_runinput_python = CODE_PATH + '/PythonSubmissionScripts/runPythonInput.sh'

    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""

    input1= "highResTimeTraceAnalysisTriangle4.py 4000 " + fullPath
    input2= "multipleAVIFlash('"+ fullPath +"')"
    qsubCommand1 = ("sbatch --mem=4000 " 
        + "--time=3:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/datFlash-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/datFlash-%J.err" + "\""
        + " " + code_runinput_python + " " 
        + " \"" + input1 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=2000 "  
        + "--time=2:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/avFlash-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/avFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input2 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    return commandList
    
def flashNewNotWorking_input(commandList,fullPath, email_flag = False):
    commandList.insert(len(commandList)-1, '####TIME SYNC####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'
    
    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""
        
    memlimit = 4000

    input1= "python highResTimeTraceAnalysisTriangle4.py 4000 " + fullPath
    input2= "multipleAVIFlash('"+ fullPath +"')"
    qsubCommand1 = ("sbatch --mem=4000 " 
        + "--time=3:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/datFlash-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/datFlash-%J.err" + "\""
        + input1 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=2000 "  
        + "--time=2:02:00" 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/avFlash-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/avFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input2 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    return commandList
    

def custom_input(commandList,input_command,fullPath, email_flag = False,time='180',mem='16000'):
    commandList.insert(len(commandList)-1, '####CUSTOM INPUT####'+NOW)
    folderName=get_folder_name(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    time_str=" --time=" + time
    mem_str=" --mem=" + mem
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runCustomMatlabInput.sh'
    
    if email_flag:
        email_script=get_email_script('end,fail')
    else:
        email_script=""
        
    qsubCommand = ("sbatch"
        + mem_str 
        + time_str 
        + " -D " + folderName
        + " -J "+ folderName 
        + email_script
        + " --output=\"" + outputFilePath + "/custom-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/custom-%J.err" + "\""
        + " " + code_runinput 
        + " \"" + input_command +"\" ")
    
    commandList.insert(len(commandList)-1, qsubCommand)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


#Write all of the inputs that were submitted into della into a text file and place it in the output folder. This file can be copied directly into the terminal of della to re run the job. 
def write_input(commandList,client,fullPath):
    
    #also adding write comments here, could be almost anywhere 
    write_comments(fullPath)
    
    outputFilePath=make_output_path(fullPath)
    fileName=outputFilePath+'/input.txt'
#open sftp client to do the write, this is needed for writing from local machine over ssh, otherwise, just write normally. 
    if socket.gethostname()=='tigressdata.princeton.edu' or socket.gethostname()=='tigressdata2.princeton.edu':
        print(fileName)
        with open(fileName,'a') as f:
            for command in commandList:
                f.write(command)
                f.write('\r\n')
            os.chmod(fileName,06775)
    else:
        ftp = client.open_sftp()
        try:
            ftp.mkdir(outputFilePath)
        except IOError:
            pass
        file=ftp.open(fileName, 'a')
        for command in commandList:
            file.write(command)
            file.write('\r\n')
        file.flush()
        ftp.close()
        
        
#search for the comments file in the parent directory and parse the comments pertaining to the folder being analyzed. Save that output into the folder being analyzed.

def write_comments(fullPath):
    
    #only works on tigress for now
    if not socket.gethostname()=='tigressdata.princeton.edu':
        print('Not on tigress! comments file not writing')
        return
    
    #don't overwrite if data_comments file already exists
    data_comments_file = fullPath + '/data_comments.txt'
    if os.path.exists(data_comments_file):
        print('Comments file already exists! not overwriting')
        return


    parent_dir,folder_name = os.path.split(fullPath)
    comments_file=parent_dir+'/Comments.txt'

    
    if os.path.exists(comments_file):
        
        target=False
        w=open(data_comments_file,'w')

        
        with open(comments_file,'r') as f:
            search_line=f.readlines()
            for line in search_line:
            #write lines between seeing the date and seeing the word 'NEW FILE'
                if folder_name in line:
                    target=True
                if 'NEW FILE' in line and target:
                     w.close()
                     break
                if target:
                     w.write(line)
    else: 
        print('Parent folder comments file not found')


#make the outputfolder where things are going to live, normally this folder is created automatically by the jobs but we need to make it ourselves because thats where we chose to put the input file from write_input. 
def make_ouputfolder(client,fullPath):
    outputFilePath=make_output_path(fullPath)
    
    if socket.gethostname()=='tigressdata.princeton.edu':
        if not os.path.exists(outputFilePath):
            os.makedirs(outputFilePath)
            os.chmod(outputFilePath,06775)
    else:
        ftp = client.open_sftp()
        try:
            ftp.chdir(outputFilePath) # sub-directory exists
        except IOError:
            ftp.mkdir(outputFilePath)
            time.sleep(1)
            ftp.chmod(outputFilePath,06775)
            ftp.close()
