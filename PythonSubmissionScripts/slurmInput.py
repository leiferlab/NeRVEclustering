#!/usr/bin/python
# Supplementary code to go with submitWormAnalysis code developed in the Leifer Lab. 

# Jeffrey Nguyen

import os
import numpy as np
import datetime
import pickle
import subprocess
import getpass

CODE_PATH='/tigress/LEIFER/communalCode/3dbrain'
qString_min = "--time=180"
PS_NAME1 =  'PointsStats.mat'
PS_NAME2 =  'PointsStats2.mat'
NOW=datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y")

def get_git_hash(commandList,client):
    git_command='cd '+ CODE_PATH +'\n git rev-parse HEAD \n cd $HOME'
    stdin,stdout,stderr = client.exec_command(git_command)
    current_git_hash=stdout.readlines()
    commandList.insert(len(commandList)-1, '#### current git: '+current_git_hash[0])
    return commandList
    
def path_setup(commandList):
    commandList.insert(len(commandList)-1, '####PATH SETUP####'+NOW)
    code_home,_=os.path.split(CODE_PATH)
    commandList.insert(len(commandList)-1, "export CODE_HOME="+code_home)
    commandList.insert(len(commandList)-1, "umask 002")
    return commandList

def pickle_load():
    # get ready for pickled variables 
    pickle_path = (os.path.expanduser('~') + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    if not os.path.exists(pickle_path):
            os.makedirs(pickle_path)
        
    if not os.path.exists(pickle_file):
            storedUsername = { "username": "USER" }
            pickle.dump( storedUsername, open(pickle_file, "wb" ) )
    
    # check to see if there is a username in the pickle file
    prevUser = pickle.load( open( pickle_file, "rb" ) )
    return prevUser
    
def get_email_script():
    user= getpass.getuser()
    user_email=user+'@princeton.edu'
    mail_script= ' --mail-type=end --mail-user=' + user_email
    return mail_script

def make_output_path(fullPath):
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    return outputFilePath
    

def centerline_input(commandList,fullPath,email_flag = False):
    commandList.insert(len(commandList)-1, '####CENTERLINES####'+NOW)
    
    folderName=os.path.basename(fullPath)
    outputFilePath= make_output_path(fullPath)
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    code_centerline = CODE_PATH + '/PythonSubmissionScripts/runWormCenterlineFitting.sh'
    code_centerline_compile = CODE_PATH + '/PythonSubmissionScripts/runWormCenterlineCompile.sh'
    
    input0 = "clusterStraightenStart('"+ fullPath + "')"
    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""
    
    qsubCommand0 = ("sbatch --mem=12000 " 
        + qString_min 
        + " -D " + folderName
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/CLstart.out"+"\""
        + " --error=\"" + outputFilePath + "/CLstart.err" + "\""
        + " " + code_runinput
        + " \"" + input1 +"\" ")
        
    commandList.insert(len(commandList)-1, qsubCommand0)
    
    qsubCommand1 = ("sbatch --mem=12000 " 
        + qString_min 
        + " -D " + folderName
        + " -J "+ folderName
        + " -d singleton"
        + " --output=\"" + outputFilePath + "/CLjob-%J.out"+"\""
        + " --error=\"" + outputFilePath + "/CLjob-%J.err" + "\""
        + " --array=1-32:1"
        + " " + code_centerline
        + " '"  + fullPath +"' ")
        
    commandList.insert(len(commandList)-1, qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=12000 " 
        + qString_min 
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
    folderName=os.path.basename(fullPath)
    outputFilePath= make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'
    code_straighten = CODE_PATH + '/PythonSubmissionScripts/runWormStraighten.sh'
    code_pscompiler = CODE_PATH + '/PythonSubmissionScripts/runWormCompilePointStats.sh'
    
    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""
    
    nRuns=np.ceil(totalRuns/1000)
    stepSize=totalRuns//300
    
    input0 = "clusterStraightenStart('"+ fullPath + "')"
    qsubCommand0 = ("sbatch --mem=12000 " 
        + qString_min 
        + " -D " + folderName
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/straight_s-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/straight_s-%J.err" + "\""
        + " " + code_runinput
        +" \"" + input0 +"\"")
    qsubCommand0 = "q0=$("+ qsubCommand0 + ")"
    commandList.insert(len(commandList)-1, qsubCommand0)
    commandList.insert(len(commandList)-1, "echo $q0")
    commandList.insert(len(commandList)-1, '\r')

    dependencyString=" --dependency=afterok:${q0##* }"

    
    for i in range(0,int(nRuns)):
        offset=str(i*1000)
        if i>=(nRuns-1):
            currentLimit=str(int(totalRuns)%1000)
        else:
            currentLimit="1000"

        qsubCommand1 = ("sbatch --mem=12000 " 
            + qString_min + " -D " + folderName
            + " -J "+ folderName 
            + dependencyString
            + " --output=\"" + outputFilePath + "/straight-%J.out" 
            + "\" --error=\"" + outputFilePath + "/straight-%J.err" + "\""
            + " --array=1-" + currentLimit + ":" + str(stepSize) 
            + " " + code_straighten 
            + " '"  + fullPath +"' "  + str(stepSize) +" " + offset)
        commandList.insert(len(commandList)-1, qsubCommand1)
    commandList.insert(len(commandList)-1, '\r')
        
    qsubCommand2 = ("sbatch --mem=12000 " 
        + qString_min 
        + " -D " + folderName
        + " -J "+ folderName 
        + " -d singleton"
        + email_script
        + " --output=\"" + outputFilePath + "/pscompile-%J.out" + "\" "
        + "--error=\"" + outputFilePath + "/pscompile-%J.err" + "\""
        + " " + code_pscompiler 
        +" '" + fullPath +"'")
    commandList.insert(len(commandList)-1, qsubCommand2)
    commandList.insert(len(commandList)-1, '\r')
    return commandList



def track_input(commandList,fullPath,totalRuns,nRef,email_flag = False):
    commandList.insert(len(commandList)-1, '####TRACKING####'+NOW)
    totalRuns=int(totalRuns)
    nRef=int(nRef)
    
    folderName=os.path.basename(fullPath)
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    
    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    code_track = CODE_PATH + '/PythonSubmissionScripts/runWormCellTracking.sh'
    code_trackcompiler = CODE_PATH+'/PythonSubmissionScripts/runWormTrackCompiler.sh'
    
    
    qString_track = "--time=" + str(np.max((nRef*2,180)))     

    matlabDirName = fullPath + "/" +  PS_NAME1
    matlabDirName2 = fullPath + "/" + PS_NAME2
    
    stepSize=int(np.ceil(50.0/nRef))
    
    input1= "makePointStatsRef('"+ fullPath +"',"+ str(nRef) + ")"
    qsubCommand0 = ("sbatch --mem=2000 " 
        + qString_min 
        + " -D " + folderName
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/track_s-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/track_s-%J.err" + "\""
        + " " + code_runinput
        +" \"" + input1 +"\"")
    qsubCommand0 = "q1=$("+ qsubCommand0 + ")"
    commandList.insert(len(commandList)-1, qsubCommand0)
    commandList.insert(len(commandList)-1, "echo $q0")
    commandList.insert(len(commandList)-1, '\r')
    dependencyString=" --dependency=afterok:${q1##* }"
    
    nRuns=totalRuns//1000+1

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
    
    qsubCommand2 = ("sbatch --mem=100000 " 
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
    
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    if email_flag:
        email_script=get_email_script()
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
        
    qsubCommand6 = ("sbatch --mem=100000 " 
        + qString_min 
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
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    
    input1= "fiducialCropper3('"+ fullPath +"')"

    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""
        
    qsubCommand7 = ("sbatch --mem=16000 " 
        + qString_min 
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
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'

    if email_flag:
        email_script=get_email_script()
    else:
        email_script=""

    input1= "highResTimeTraceAnalysisTriangle4('"+ fullPath + "')"
    input2= "multipleAVIFlash('"+ fullPath +"')"
    qsubCommand1 = ("sbatch --mem=2000 " 
        + qString_min 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/datFlash-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/datFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input1 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=2000 "  
        + qString_min 
        + " -J "+ folderName
        + email_script
        + " --output=\"" + outputFilePath + "/avFlash-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/avFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input2 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand2)
    print(qsubCommand2)
    return commandList
    
    
def write_input(commandList,client,fullPath):
    outputFilePath=make_output_path(fullPath)
    fileName=outputFilePath+'/input.txt'
    
    ftp = client.open_sftp()
    file=ftp.file(fileName, "a", -1)
    for command in commandList:
        file.write(command)
        file.write('\r\n')
    file.flush()
    ftp.close()
    
def make_ouputfolder(client,fullPath):
    outputFilePath=make_output_path(fullPath)
    ftp = client.open_sftp()
    try:
        ftp.chdir(outputFilePath) # sub-directory exists
    except IOError:
        ftp.mkdir(outputFilePath)
    ftp.close()
    