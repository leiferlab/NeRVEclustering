#!/bin/sh
#
# File:   wormAnalysis_makefile.sh
# Author: jnguyen

# This program sets up paths and ssh keys for running code from the 3dbrain repo. 
# this is made to run from terminal in tigressdata. This be done by ssh or VNC
# For SSH, ssh into USER@tigressdata.princeton.edu.

# For VNC, follow the instructions from https://www.princeton.edu/researchcomputing/faq/how-do-i-use-vnc-on-tigre/
# to set up a VNC window. Open the terminal and run this script



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
	
	pass2=$(ssh $USER@della.princeton.edu -qo PasswordAuthentication=no echo 0 || echo 1)

	if [ "$pass2" == "1" ]; then
		echo "ERROR: still have problems connecting to della without password"
	else
		echo "Success! Keys saved for della connection"
	fi
	
	# load python module, install paramiko with pip
	module load anaconda
	pip install --user paramiko
	
else
	echo "not currently on tigressdata, make will not run"
fi

