 
Worm analysis protocol:


This repository hold the code used for the analyzing movies from the Leifer Lab's Whole brain imaging set. The details of the pipeline are described in the paper “Automatically tracking neurons in a moving and deforming brain” by Nguyen et al 2017.  The paper can be found at https://doi.org/10.1371/journal.pcbi.1005517



#########################################################################
SETUP

#########################################################################

All of the analysis is done in matlab, but many of them are called on DELLA, which is Princeton University’s SLURM based computational cluster. Jobs are submitted to della via python wrappers that take in some inputs. Folders with HighMag data are on tigress. The corresponding low mag folder should be placed inside the high mag folder. Prior to running submission scripts, you need to have access to della,tigressdata, /tigress/LEIFER (ask Andy to email John Wiggins). If you are using a Windows machine, you will need to download and install PUTTY. 



TO RUN THE CODE FROM TIGRESSDATA


Once you have access to della and tigressdata, open a VNC connection by following the instructions on https://www.princeton.edu/researchcomputing/faq/how-do-i-use-vnc-on-tigre/. 

Open a terminal window by going to Applications->System Tools -> Terminal

Run the following commands:



/tigress/LEIFER/communalCode/3dbrain/PythonSubmissionScripts/wormAnalysis_makefile.sh 
	


#####TO RUN A PYTHON SCRIPT#####:



Most of the anlaysis is run via Python submission scripts. From Terminal navigate to the code location with:

	

cd /tigress/LEIFER/communalCode/3dbrain/PythonSubmission/



You can then run the Python submission codes by entering:

	module load anaconda

python <python submission code name>.py
	



#####TO RUN A MATLAB SCRIPT#####:


If running from tigressdata, matlab can be found by typing this into terminal:

	/usr/licensed/matlab-R2017a/bin/matlab



In the matlab command line, set up the paths to use these programs with:



	cd /tigress/LEIFER/communalCode
	
	path(pathdef)



You can then run any GUI typing the name of the .m file into the command line.



######################################################################

ANALYSIS PIPELINE

#####################################################################

To start, you must put the LowMag folder which has the .avi videos inside the corresponding BrainScanner folder. Run each of the Python submission codes as described above. 


STEP 0a: TIMING SYNCHRONIZATION FOR VIDEOS

	Python submission code:

		submitWormFlashFinder.py

	Matlab analysis code called by python:

		highResTimeTraceAnalysisTriangle4.m

		multipleAVIFlash.m


	File Outputs:

			hiResData.mat file -data for timing of .dat file

			cam0flashtrack.mat

			cam1flashtrack.mat



STEP 0b: IMAGE ALIGNMENT FOR VIDEOS
(done on tigressdata VNC with Matlab, does not depend on timing)

After taking the alignment videos on both computers, move the LowMag folder into the the BrainScanner folder. This is likely done on the computer "Bardeen" or on tigressdata VNC. Use alignment_gui.m on the BrainScanner folder that has the alignment videos. After saving the alignments, move the alignment.mat file into each of the BrainScanner folders for analysis. 


STEP 1: WORM CENTERLINE DETECTION
(done locally or on tigressdata VNC for manual centerline initialization)

	wormCL_tip_clicker.m
		- this is an optional GUI that will allow the user to help the centerline fitting by explicitly clicking on the location of the head and the tail. 
	initializeCLWorkspace.m 
		- this will calculate background images and allow the user to manually initialize a few centerlines to help the detection algorithm.
 
		- when this is done, move the alignment and the initialCLWorkspace.mat into the corresponding LowMag folder on tigress.

 
	Python submission code:

		submitWormAnalysisCenterline.py
	Matlab analysis code:

		clusterWormCenterline.m
	File Outputs:
	CLstartworkspace.mat, initialized points and background images for darkfield images
	CL_files folder, containing partial CL.mat files
	BehaviorAnalysis folder, containing the centerline.mat file with XY coordinates for each image.


	*NOTE: due to poor image quality of dark field images, it may be necessary to use some of the code developed by ANL to manually adjust centerlines



**steps 2-5 all use submitWormAnalysisPipelineFull.py for submission. 
STEP 2: STRAIGHTEN AND SEGMENTATION
	
Python submission code:

		submitWormAnalysisPipelineFull.py


	Matlab analysis code called by python:


		clusterStraightenStart.m

		clusterWormStraightening.m
	
File Outputs:
	startWorkspace.mat, initial workspace used for during straightening for all volumes
	CLStraight* folder, folder containing all saved straightened tif files and results of segmentation.


STEP 3: NEURON REGISTRATION VECTOR ENCODING AND CLUSTERING

	Python submission code:

		submitWormAnalysisPipelineFull.py


	Matlab analysis code called by python:



		clusterWormTracker.m

		clusterWormTrackCompiler.m

	File Outputs:
		TrackMatrixFolder, containing all registrations of sample volumes with reference volumes.
		pointStats.mat, structure containing all coordinates from all straightened volumes along with a trackIdx, the result of initial tracking of points. 

STEP 4: ERROR CORRECTION

	Python submission code:

		submitWormAnalysisPipelineFull.py


	Matlab analysis code called by python:

		clusterBotChecker.m

		clusterBotCheckCompiler.m

	File Outputs:
 	botCheckFolder- folder containing all coordinate guesses for all times, one mat file for each neuron.
	pointStatsNew.mat- matfile containing the refined trackIdx after error correction. 

STEP 5: SIGNAL EXTRACTION

	Python submission code:

		submitWormAnalysisPipelineFull.py


	Matlab analysis code called by python:


		fiducialCropper3.m

	File Output:
	heatData.mat, all signal results from extracting signal from the coordinates. 
	

#########################################################################
USEFUL GUIS FOR VISUALIZATION


These guis use wasd gaming controls, ie a for left, d for right, w for up, s for down:

GUI to Check Data 
Before any analysis

------------------------------------
ScanBinaryImageStack.m - Gui to view raw .dat file movies. This also works with .avi files



GUIs for post-centerlines

-------------------------
wormCLviewer.m - Gui to view darkfield worm images along with the centerline


WormAnalysisPreview.m - Gui to check time and spatial alignments of all videos. Good for use post centerlines but prior to straightening.


GUIs for post-pipeline to make sure things worked
-------------------------------------------------

VisualzeWorm3danalysis.m - Check that behavior and straightening worked well. Also shows tracked neurons. 

VisualizeTrackedData.m - Check to see that tracking works well. Works on unstraightened .dat file. 





#########################################################################
FOLDER DATA REQUIREMENTS



======Worm Videos======


sCMOS_Frames_U16_1024x1024.dat  -	binary image file from the Leifer Lab's WholeBrainImaging labview program. The images do not need to be 1024 by 1024.The image is produced by a dual view setup so that the left and right images represent the RFP and the GCaMP6s images in some order. 


LowMagBrain* folder containing all low magnification data including

cam0.avi	-	low magnification fluorescent images of the worm’s brain

cam1.avi	-	low magnification dark field images of the worm’s posture

CamData.txt	-	text file with relative timing for every frame



*****NOTE: FOR OLDER DATA (pre 2016)
We used to use avi’s that were not time synced and used a YAML file containing the meta data. Using this data requires the code from the repo https://github.com/leiferlab/MindControlAccessUtils.git.
****** (end note for older data)
 

======Raw Text files======


These files contain timing information for every frame of each of the video feeds. They also contain information about the positions of the stage and the objective. This information, along with the videos themselves, are used to align the timing for all of the videos. Several camera flashes are used throughout the recording. The timing of  

labJackData.txt - 	Raw outputs from LabVIEW for the stage, the piezo that drives the objective, the sCMOS camera, and the function generator (FG), taken at 1kHz. The objective is mounted on a piezo stage that is driven by the output voltage of the function generator. The 1kHz clock acts as the timing for each event.

	Columns:
	FuncGen - 	Programmed output from FG, a triangle wave at 6 Hz. 
	Voltage - 	Actual FG output
	Z Sensor - 	voltage from piezo, which controls Z position of objective.
	FxnGen Sync- 	Trigger output from FG, triggers at the center of the triangle.
	Camera Trigger-	Voltage from HiMag camera, down sweeps indicate a frame has been grabbed from HiMag Camera
	Frame Count - 	Number of frames that have been grabbed from HiMag Camera. Not all frames that are grabbed are saved, and the saved frames will be 			indicated in the saved frames field in the next text file. 
	Stage X -	X position from stage
	Stage Y - 	Y position from stage


CameraFrameData.txt - 	Metadata from each grabbed frame for HiMag  images, saved in LabVIEW. The timing for each of these frames can be pulled from the labJackData.txt.

	Columns:
	Total Frames -	Total number of grabbed frames
	Saved Frames - 	The current save index, not all grabbed frames are saved. If this increments, the frame has been saved.
	DC offset    -	This is the signal sent to the FG to translate center of the triangle wave to keep the center of the worm in the middle of the wave.
	Image STdev - 	standard deviation of all pixel intensities in each image.



alignments.mat -matlab structure with fields containing the affine transforms for image alignment. Alignments created using the createAlignment.m program: The file has the following fields
	lowResFluor2BF	- affine transform between low magnification fluorescent and behavior images
	S2AHiRes	- affine transformation between the RFP channel and the Gcamp6s channel of the hi magnification image along the the cropping rectangles for each of them.
	Hi2LowResFluor	- affine transformation between the RFP hi magnification image (before cropping) and the hi magnification fluorescent images



======Output Files ========
The result of the pipeline is various .mat files that contain neural signals and coordinates. The main output is heatData.mat

heatData.mat 
Signal vairables
	rRaw - an N neurons x T volumes matrix with the raw red signal from each of the neurons. Signals are averaged pixel values around each tracked neuron with no other processing except for flash removal. 
	gRaw - same as rRaw but for the green signal. 

	rPhotoCorr - the rRaw signal after photobleaching correction for each neuron. No other smoothing or normalization is applied. Photobleaching correction is applied by fitting an exponential curve to a wide 20th percentile filter, and then subtraction the exponential from the raw signal.
	gPhotoCorr - same as above but with the green signal. Exponential curves are fit independently. 

	R2 - Smoothed and normalized version of rPhotoCorr.Normalization is done as delta F/ F0, where F0 is the lower 20th percentile signal. A 5 time step (.83s) fwhm Gaussian is used to smooth the result. 
	G2 - Same as above but with gPhotoCorr.

	Ratio2 - First, both the green and red signals are smoothed with a 5 time step fwhm Gaussian filter. Then the Ratio is then taken as  gPhotoCorr/rPhotoCorr and normalized as delta R/ R0 in the same way as R2 and G2. Some modifications are made to deal with some of the nans. 



Other fields:
behavior - structure with 
    ethogram: t Volumes by 1 vector of behaviors, -1 for reverse, 0 pause, 1 forward, 2 turn. Behaviors determined automatically using the centerlines.
       x_pos: t Volumes by 1 vector of x coordinates in the reference frame of the plate.
       y_pos: t Volumes by 1 vector of y coordinates in the reference frame of the plate.
           v: t Volumes by 1 vector of worm center of mass velocities in the reference frame of the plate. Positive is forward, negative reverse.
       pc1_2: t Volumes by 2 vector of the projections onto the first two eigenworms. 
        pc_3: t Volumes by 1 vector of the projections onto the third eigenworm. 

XYZcoord - N neurons x 3 XYZ coordinates for the neurons. The coordinates are taken by a random straightened volume from the recording to serve as an example.

acorr - a n Neurons x n Neurons pearson correlation matrix of Ratio2

cgIdx - ordered indices derived from heirarchically clustering the correlation matrix. To show organized traces, use Ratio2(cgIdx,:)

hasPointsTime- a t volumes by 1 vector of time for each volume. 


