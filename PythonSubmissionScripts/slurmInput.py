#!/usr/bin/python
"""

"""
import os
import numpy as np
import datetime

CODE_PATH='/tigress/LEIFER/communalCode/3dbrain'
qString_min = "--time=180"
PS_NAME1 =  'PointsStats.mat'
PS_NAME2 =  'PointsStats2.mat'


def path_setup(commandList):
    code_home,_=os.path.split(CODE_PATH)
    commandList.insert(len(commandList)-1, "export CODE_HOME="+code_home)
    return commandList

def make_output_path(fullPath):
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate
    return outputFilePath


def straighten_input(commandList,fullPath,totalRuns):
    commandList.insert(len(commandList)-1, '####STRAIGHTENING####')

    totalRuns = int(totalRuns)
    folderName=os.path.basename(fullPath)
    outputFilePath= make_output_path(fullPath)
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    code_straighten = CODE_PATH + '/PythonSubmissionScripts/runWormStraighten.sh'
    code_pscompiler = CODE_PATH + '/PythonSubmissionScripts/runWormCompilePointStats.sh'
    
    nRuns=np.ceil(totalRuns/1000)
    stepSize=totalRuns//300
    
    input0 = "clusterStraightenStart('"+ fullPath + "')"
    qsubCommand0 = ("sbatch --mem=2000 " 
        + qString_min + " -D " + folderName
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

        qsubCommand1 = ("sbatch --mem=2000 " 
            + qString_min + " -D " + folderName
            + " -J "+ folderName 
            + dependencyString
            + " --output=\"" + outputFilePath + "/straight-%J.out" 
            + "\" --error=\"" + outputFilePath + "/straight-%J.err" + "\""
            + " --array=1-" + currentLimit + ":" + str(stepSize) 
            + " " +code_straighten 
            + " '"  + fullPath +"' "  + str(stepSize) +" " + offset)
        commandList.insert(len(commandList)-1, qsubCommand1)
    commandList.insert(len(commandList)-1, '\r')

        
    qsubCommand2 = ("sbatch --mem=2000 " 
        + qString_min + " -D " + folderName
        + " -J "+ folderName + " -d singleton"
        + " --output=\"" + outputFilePath + "/pscompile-%J.out" + "\" "
        + "--error=\"" + outputFilePath + "/pscompile-%J.err" + "\""
        + code_pscompiler 
        +" '" + fullPath +"'")
    commandList.insert(len(commandList)-1, qsubCommand2)
    commandList.insert(len(commandList)-1, '\r')
    return commandList



def track_input(commandList,fullPath,totalRuns,nRef):
    commandList.insert(len(commandList)-1, '####TRACKING####')
    totalRuns=int(totalRuns)
    nRef=int(nRef)
    
    folderName=os.path.basename(fullPath)
    outputFilePath= fullPath + "/outputFiles"
    currentDate=datetime.date.today()
    currentDate=str(currentDate)
    outputFilePath= outputFilePath + currentDate

    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    code_track = CODE_PATH + '/PythonSubmissionScripts/runWormCellTracking.sh'
    code_trackcompiler = CODE_PATH+'/PythonSubmissionScripts/runWormTrackCompiler.sh'
    
    
    qString_track = "--time=" + str(np.max((nRef*2,180)))     

    matlabDirName = fullPath + "/" +  PS_NAME1
    matlabDirName2 = fullPath + "/" + PS_NAME2
    
    
    
    input1= "makePointStatsRef('"+ fullPath +"',"+ str(nRef) + ")"
    qsubCommand0 = ("sbatch --mem=2000 " 
        + qString_min + " -D " + folderName
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/straight_s-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/straight_s-%J.err" + "\""
        + " " + code_runinput
        +" \"" + input1 +"\"")
    qsubCommand0 = "q0=$("+ qsubCommand0 + ")"
    commandList.insert(len(commandList)-1, qsubCommand0)
    commandList.insert(len(commandList)-1, "echo $q0")
    commandList.insert(len(commandList)-1, '\r')
    dependencyString=" --dependency=afterok:${q0##* }"
    

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
            + " --array=1-" + currentLimit
            + " " + code_track + " '" 
            + fullPath +"' " + offset)
        
        commandList.insert(len(commandList)-1, qsubCommand1)
    commandList.insert(len(commandList)-1, '\r')
    
    qsubCommand2 = ("sbatch --mem=100000 " 
        + qString_track + " -D " + folderName
        + " -J "+ folderName + " -d singleton"
        + " --output=\"" + outputFilePath + "/trackCompiler-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/trackCompiler-%J.err" + "\""
        + " " + code_trackcompiler 
        +" '" + matlabDirName +"' '" + matlabDirName2 + " '")
    commandList.insert(len(commandList)-1, qsubCommand2)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


def check_input(commandList,fullPath,nCheck,nNeurons):
    commandList.insert(len(commandList)-1, '####CHECKING####')
    nCheck=int(nCheck)
    nNeurons=int(nNeurons)
    
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_checkcompiler = CODE_PATH + '/PythonSubmissionScripts/runWormBotCheckCompiler.sh'
    code_check = CODE_PATH+'/PythonSubmissionScripts/runWormBotChecker.sh'
    qString_check = "--time="+ str(np.max((nCheck,180)))   

    matlabDirName2 = fullPath + "/" + PS_NAME2
    
    qsubCommand5 = ("sbatch --mem=8000 " 
        + qString_check + " -D " + folderName
        + " -J "+ folderName + " -d singleton"
        + " --output=\"" + outputFilePath + "/check-%J.out" + "\" "
        + " --error=\"" + outputFilePath + "/check-%J.err" + "\" "
        + " --array=1-"+ str(nNeurons) + ":" + "1"
        + " " + code_check + " '" 
        + matlabDirName2 +"' 1  1")
    commandList.insert(len(commandList)-1, qsubCommand5)
        
    qsubCommand6 = ("sbatch --mem=100000 " 
        + qString_min + " -D " + folderName
        + " -J "+ folderName + " -d singleton"
        + " --output=\"" + outputFilePath + "/checkcompile-%J.out" + " \""
        + " --error=\"" + outputFilePath + "/checkcompile-%J.err" + " \""
 #       + " --mail-type=end" + " --mail-user=" + userEmail
        + " " + code_checkcompiler 
        + "  '" + fullPath +"' ")
    commandList.insert(len(commandList)-1, qsubCommand6)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


def crop_input(commandList,fullPath):
    commandList.insert(len(commandList)-1, '####CROPPING####')
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH+ '/PythonSubmissionScripts/runMatlabInput.sh'
    
    input1= "fiducialCropper3('"+ fullPath +"')"
    
    qsubCommand7 = ("sbatch --mem=16000 " 
        + qString_min + " -D " + folderName
        + " -J "+ folderName + " -d singleton"
        + " --output=\"" + outputFilePath + "/crop-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/crop-%J.err" + "\""
        + " " + code_runinput 
        + " \"" + input1 +"\" ")
    
    commandList.insert(len(commandList)-1, qsubCommand7)
    commandList.insert(len(commandList)-1, '\r')
    return commandList


def flash_input(commandList,fullPath):
    commandList.insert(len(commandList)-1, '####TIME SYNC####')
    folderName=os.path.basename(fullPath)
    outputFilePath=make_output_path(fullPath)
    
    code_runinput = CODE_PATH + '/PythonSubmissionScripts/runMatlabInput.sh'
    
    input1= "highResTimeTraceAnalysisTriangle4('"+ fullPath + "')"
    input2= "multipleAVIFlash('"+ fullPath +"')"
    
    qsubCommand1 = ("sbatch --mem=2000 " 
        + qString_min 
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/datFlash-%J.out"+ "\"" 
        + " --error=\"" + outputFilePath + "/datFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input1 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand1)
    print(qsubCommand1)
    
    qsubCommand2 = ("sbatch --mem=2000 "  
        + qString_min 
        + " -J "+ folderName
        + " --output=\"" + outputFilePath + "/avFlash-%J.out" + "\""
        + " --error=\"" + outputFilePath + "/avFlash-%J.err" + "\""
        + " " + code_runinput + " " 
        + " \"" + input2 +"\" ")
    commandList.insert(len(commandList)-1, qsubCommand2)
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
    