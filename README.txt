 
Worm analysis protocol:


#########################################################################
Intro

#########################################################################


This repository hold the code used for the analyzing movies from the Leifer Lab's Whole brain imaging set. The details of the pipeline are described in the paper “Automatically tracking neurons in a moving and deforming brain” by Nguyen et al 2017.  The paper can be found at https://doi.org/10.1371/journal.pcbi.1005517. The data associated with this code can be found at the IEEE DataPort repository http://dx.doi.org/10.21227/H2901H



All of the analysis is done in matlab, but many of them are called on DELLA, which is Princeton University’s SLURM based computational cluster. For public usage, some paths and commands will need to be changed in order to run on local machines or home institution computing clusters. You can run this code on the example files provided at 

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

After taking the alignment videos on both computers, move the LowMag folder into the the BrainScanner folder. . Use alignment_gui.m on the data folder that has the alignment videos. After saving the alignments, move the alignment.mat file into each of the data folders for analysis. 


STEP 1: WORM CENTERLINE DETECTION

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

EARLY OUTPUT FILES

Along the way, the pipeline produces some useful alignment files for timing and image alignment. 

hiResData.mat - 	data for each frame of HiMag video. This contains the information about the imaging plane, position of the stage, timing, and which volume each frame belongs to. 
Fields:
	Z - 		z voltage from piezo for each frame, indicating the imaging plane. 
	frameTime -	time of each frame after flash alignment in seconds.
	stackIdx -	number of the stack each frame belongs to. Each volume recorded is given an increasing number starting at 1 for the first volume. For example, the first 40 images will belong to stackIdx=1, then the next 40 will have stackIdx=2 etc etc…
	imSTD - 	standard dev of each frame
	xpos and ypos - stage position for each frame
	flashLoc -	index of the frame of each flash

*note some of these fields have an extra point at the end, just remove it to make everything the same size


*flashTrack.mat - 1xN vector where N is the number of frames in the corresponding video. The values of flashTrack are the mean of each image. It will show a clear peak when a flash is triggered. This can be used to align the videos. 

*YAML.mat - 1xN vector where N is the number of frames in the corresponding video. Each element of the mcdf has all of the metadata for each frame of the video. Using this requires code from https://github.com/leiferlab/MindControlAccessUtils.git github repo. 




alignments.mat - 	set of affine transformations between videos feeds. Each has a "tconcord" field that works 
with matlab’s imwarp function.
Fields:
	lowresFluor2BF-	Alignment from low mag fluorescent video to low mag behavior video
	S2AHiRes -	Alignment from HiMag Red channel to HiMag green channel. This alignment is prior to cropping of the HiMag Red channel. 
	Hi2LowResF -	Alignment from HiMag Red to low mag fluorescent video
	


FINAL OUTPUT FILES

The result of the pipeline is various .mat files that contain neural signals and coordinates. The main output is heatData.mat

heatData.mat 
Signal vairables
	rRaw - an N neurons x T volumes matrix with the raw red signal from each of the neurons. Signals are averaged pixel values around each tracked neuron with no other processing except for flash removal. 
	gRaw - same as rRaw but for the green signal. 

	rPhotoCorr - the rRaw signal after photobleaching correction for each neuron. No other smoothing or normalization is applied. Photobleaching correction is applied by fitting an exponential curve to a wide 20th percentile filter, and then subtraction the exponential from the raw signal.
	gPhotoCorr - same as above but with the green signal. Exponential curves are fit independently. 

	R2 - Smoothed and normalized version of rPhotoCorr.Normalization is done as delta F/ F0, where F0 is the lower 20th percentile signal. 
	G2 - Same as above but with gPhotoCorr.

	Ratio2 - The ratio signal is defined as gPhotoCorr/rPhotoCorr, the Ratio is then normalized as delta R/ R0. is the same way as R2 and G2. 



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


