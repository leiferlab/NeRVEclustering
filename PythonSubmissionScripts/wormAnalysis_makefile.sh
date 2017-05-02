#!/bin/sh
#
# File:   wormAnalysis_makefile.sh
# Author: jnguyen

# This program sets up paths and ssh keys for running code from the 3dbrain repo. 
# this is made to run from terminal in tigressdata. This be done by ssh or VNC
# For SSH, ssh into USER@tigressdata.princeton.edu. A key file, currently located in 
#/tigress/LEIFER/.ssh, is copied into the USER's home directory on tigress. The key is then
# added to the "authorized_keys" file in the USER's .ssh folder in DELLA. All users who run
# this program will be able to use their USER tigressdata account to SSH into all other USERs
# DELLA accounts. 

# The program also installs paramiko for Anaconda. 

# For VNC, follow the instructions from https://www.princeton.edu/researchcomputing/faq/how-do-i-use-vnc-on-tigre/
# to set up a VNC window. Open the terminal by going to Application>System Tools> Terminal

# From the terminal, run this file.



#check if we currently need keys to ssh in

if [ "$HOSTNAME" == "tigressdata.princeton.edu" ]; then
	pass=$(ssh $USER@della.princeton.edu -qo PasswordAuthentication=no echo 0 || echo 1)
	if [ "$pass" == "1" ]; then
		# copy the key file into the users home directory, and chmod to proper permissions
		mkdir -p $HOME/.ssh
		cp /tigress/LEIFER/.ssh/id_rsa $HOME/.ssh/id_rsa
		chmod 700 $HOME/.ssh/id_rsa
		
		cp /tigress/LEIFER/.ssh/id_rsa $HOME/.ssh/id_rsa.pub
		chmod 700 $HOME/.ssh/id_rsa
# copy the ssh key from keyfile into .ssh folder in della. you'll need to input your password 2x.
		echo " ###Keys not found, copying keys from /tigress/LEIFER###. 
			###You will need to enter your password###"
		ssh $USER@della.princeton.edu mkdir -p .ssh
		cat /tigress/LEIFER/.ssh/id_rsa.pub | ssh $USER@della.princeton.edu 'cat >> .ssh/authorized_keys'
		cat /tigress/LEIFER/.ssh/id_rsa | ssh $USER@della.princeton.edu 'cat >> .ssh/authorized_keys'
		

	else
		echo "Keys found"
	fi
	
	echo "Changeing default save permissions to 775"
	greq -q "umask 002" .bashrc
	if [ $? -eq 1 ]; then
		echo "umask 002" >> $HOME/.bashrc
	fi
		
	pass2=$(ssh $USER@della.princeton.edu -qo PasswordAuthentication=no echo 0 || echo 1)

	if [ "$pass2" == "1" ]; then
		echo "ERROR: still have problems connecting to della without password"
	else
		echo "Success! Keys saved for della connection"
		echo "Setting Della default permissions"
		ssh jnguyen@della.princeton.edu "grep -q \"umask 002\" $HOME/.bashrc"
		if [ $? -eq 1 ]; then
			echo "umask 002" | ssh $USER@della.princeton.edu "cat >> $HOME/.bashrc"
		fi
	fi
	# load virtualgl for matlab
	module laod virtualgl
	# load python module, install paramiko with pip
	module load anaconda
	pip install --user paramiko
	
else
	echo "not currently on tigressdata, make will not run"
fi

