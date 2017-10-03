#!/usr/bin/python

# guiHelper has functions that help build the gui and connect to della. The functions are called by submitWorm  analysis codes.
#


try:
    import Tkinter as tk
    import tkFileDialog
except ImportError:
    import tkinter as tk
    import tkinter.filedialog as tkFileDialog
import paramiko
import socket
import os
import pickle


#server jobs willb e submitted to
SERVER ='della.princeton.edu'
#location of the ssh key file shared by lefierdata
KEYPATH='/tigress/LEIFER/.ssh/id_rsa'

#load the pickle files for default values to put into fields. The pickle files store your previous entries in the pickles2.p file in ~. 
def pickle_load():
    # get ready for pickled variables 
    pickle_path = (os.path.expanduser('~') + "/platypusTemp/")
    pickle_file = pickle_path + "pickles2.p"
    if not os.path.exists(pickle_path):
            os.makedirs(pickle_path)
        
    if not os.path.exists(pickle_file):
            storedUsername = { "username": "USER" }
            pickle.dump( storedUsername, open(pickle_file, "wb" ) )
    
    # return the dictionary with all previous values stored. 
    prevUser = pickle.load( open( pickle_file, "rb" ) )
    
    #fill in default values
    if 'username' not in prevUser:
            prevUser['username'] = "USER"
    if 'time' not in prevUser:
            prevUser['time'] = "180"
    if 'mem' not in prevUser:
            prevUser['mem'] = "16000"
    if 'date' not in prevUser:
            prevUser['date'] = "testing_sets"
    if 'folderName' not in prevUser:
            prevUser['folderName'] = "Brain_working_dataset"
    if 'frameNumber' not in prevUser:
            prevUser['frameNumber'] = "1000"
    if 'refNumber' not in prevUser:
            prevUser['refNumber'] = "100"
    if 'neuronNumber' not in prevUser:
            prevUser['neuronNumber'] = "150"
    if 'checkNumber' not in prevUser:
            prevUser['checkNumber'] = "500"
    if 'matlab_command' not in prevUser:
            prevUser['matlab_command']='Put your code here!'
    if 'code_path' not in prevUser:
            prevUser['code_path']='Put the path to the .path file here!'
            

    return prevUser


#read the submissionParameters.txt file in the data folder to get the number of frames. 
def get_nframes(client,fullPath):
    
    print('Getting number of frames from GUI')
    subDataFile = fullPath + '/submissionParameters.txt'
    client2 = client.open_sftp()
    sub_data={}
    # read each line, saving key-values pairs into a dictionary
    try:
        with client2.open(subDataFile) as remote_file:
            for lines in remote_file:
                lines=lines.split(" ")
                print(lines)
                sub_data[lines[0]]=lines[1]
    except Exception:
        raise NameError("Remote File not found")
        
    if sub_data.get('NFrames') != None:
        print('Getting number of frames from Remote')
        NFrames=sub_data['NFrames']
        return NFrames
    else: 
        raise NameError("Cannot connect to get number of Frames")
    

def dellaConnect(username,password=None):
    #use the username to connect to della, if password is provided, use it, otherwise, use ssh keys
    if socket.gethostname()=='tigressdata.princeton.edu':
        #if on tigressdata, use path to key file file in /tigress/LEIFEr
        key = paramiko.RSAKey.from_private_key_file(KEYPATH)
    else:
        key=None
    # connect and submit job via sbatch
    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    if not password:
        client.connect(SERVER, 22, username,pkey=key)
    else:
        try:
            client.connect(SERVER, 22, username, password)
        except Exception:
            raise Exception('Cannont connect, problem with username password combo')
    return client
    
def passwordReqCheck(username):
    # check if a password is required for a given username
    # try to see if a password is needed
    try:
        dellaConnect(username)
        isNeedsPassword = False
    except paramiko.AuthenticationException:
        isNeedsPassword = True
    return isNeedsPassword
    
def selectFolder(master=None):
    # select a folder using a popup window
    folder=tkFileDialog.askdirectory(mustexist=False , initialdir= '/tigress/LEIFER/PanNeuronal/')
    if folder:
        #parse the folder path and populate the filed on the gui. 
        path,folderName=os.path.split(folder)
        path,date=os.path.split(path)
        master.e['parent_path'].delete(0,tk.END)
        master.e['parent_path'].insert(0,path)
        master.e['date'].delete(0,tk.END)
        master.e['date'].insert(0,date)
        master.e['folder_name'].delete(0,tk.END)
        master.e['folder_name'].insert(0,folderName)
        print(folder)
        
def selectPath(master=None):
    filename = tkFileDialog.askopenfilename(initialdir= '/tigress/LEIFER/communalCode/ ')
    master.e['code_path'].delete(0,tk.END)
    master.e['code_path'].insert(0,filename)
# 
# #class for building the gui and populating the rows
class submitTK(tk.Tk):
    #build the gui with some number of max rows and cols, 
    #submitTK is a subclass of tkinter's Tk(). 
    def __init__(self, rows=10, cols=2):
        
        tk.Tk.__init__(self)
        #build rows and cols
        for i in range(rows):
            self.columnconfigure(i,pad=3)
        for j in range(cols):
            self.rowconfigure(j, pad=3)
        # make counter for the row being populated
        self.row_count=0
        # self.e will have all the values of the rows stored
        self.e=dict()
        
        self.title("Options")
        
    # add a text row to the field, the row will have a text label, 
    # row is added to the next open row. 
    # field name for using in the e, 
    # default field entry
    # show, for passwords, use '*'
    
    def addGuiField(self,label,name,default="",show=None):
        master_label = tk.Label(self, text=label)
        master_label.grid(row=self.row_count, column=0,sticky=tk.W+tk.E)
        self.e[name] = tk.Entry(self,show=show)
        self.e[name].insert(0, default)
        self.e[name].grid(row=self.row_count, column=1,sticky=tk.W+tk.E)
        self.row_count+=1
        
    # make checkbox, again with label and field name
    def addGuiCheck(self,label,name,default=1):
        master_label = tk.Label(self, text=label)
        master_label.grid(row=self.row_count, column=0, sticky=tk.W+tk.E)
        var1= tk.IntVar()
        self.e[name]= tk.Checkbutton(self, text=None, variable=var1)
        self.e[name].var = var1
        self.e[name].grid(row=self.row_count, column=1,sticky=tk.W+tk.E)
        self.e[name].var.set(1)
        self.row_count+=1
        
    # add a button with a text label and a function callback handle. 
    def addGuiButton(self,label,b_command=None):
        self.b = tk.Button(self, text=label, width=10, command=b_command)
        self.b.grid(row=self.row_count,columnspan=2,sticky=tk.W+tk.E)
        self.row_count+=1
        
    #save values fo fields for defaults during next call
    def pickleDump(self):
        pickle_path = (os.path.expanduser('~') + "/platypusTemp/")
        pickle_file = pickle_path + "pickles2.p"
        prevUser=pickle_load()
        #refill prevUser dict with master entries
        if 'user_name' in self.e:
            prevUser['username']=self.e['user_name'].get()
        if 'time' in self.e:
            prevUser['time']=self.e['time'].get()
        if 'mem' in self.e:
            prevUser['mem']=self.e['mem'].get()
        if 'date' in self.e:
            prevUser['date'] = self.e['date'].get()
        if 'folder_name' in self.e:
            prevUser['folderName']=self.e['folder_name'].get()
        if 'nframes' in self.e:
            prevUser['frameNumber']=self.e['nframes'].get()
        if 'n_ref' in self.e:
            prevUser['refNumber']=self.e['n_ref'].get()
        if 'n_neurons' in self.e:
            prevUser['neuronNumber']=self.e['n_neurons'].get()
        if 'n_checks' in self.e:
            prevUser['checkNumber']=self.e['n_checks'].get()
        if 'matlab_command' in self.e:
            prevUser['matlab_command']=self.e['matlab_command'].get()
        if 'code_path' in self.e:
            prevUser['code_path']=self.e['code_path'].get()
            
        #save prevUser as pickle
        pickle.dump(prevUser, open(pickle_file, "wb" ) )
        
        
    # run the main loop. 
    def run(self):
        tk.mainloop()
        
