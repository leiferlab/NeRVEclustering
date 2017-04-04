 
Worm analysis protocol:


This repository hold the code used for the analyzing movies from the Leifer Lab's Whole brain imaging set. The details of the pipeline are described in the paper “Automatically tracking neurons in a moving and deforming brain” by Nguyen et al 2016. 

#########################################################################
QUICK SUMMARY

All of the analysis is done in matlab, but many of them are called on DELLA, which is Princeton University’s SLURM based computational cluster. Jobs are submitted to della via python wrappers that take in some inputs. Folders with HighMag data are on tigress. The corresponding low mag folder should be placed inside the high mag folder. Prior to running submission scripts, you need to have access to della, /tigress/LEIFER (ask Andy to email John Wiggins), save your ssh keys (see http://www.linuxproblem.org/art_9.html for mac), have python installed with paramiko.  

Inside you’re home directory on della, you need to make directories “scripts” and “data”. Inside data, you can clone the git repo found at https://github.com/leiferlab/3dbrain.git. For now, you also need shae_pythonSubmissionScripts. This is on sheavitzdata. Some of the paths are hard coded into the python scripts. 





STEP 0a: TIMING SYNCHRONIZATION FOR VIDEOS

	Python submission code:

		submitWormFlashFinder.py

	Matlab analysis code:

		highResTimeTraceAnalysisTriangle4.m

		multipleAVIFlash.m


	File Outputs: 	*YAML.mat files
			*flashTrack.mat files
			hiResData.mat file




STEP 0b: IMAGE ALIGNMENT FOR VIDEOS
(done locally for point matching)
	
Use ExtractAlignmentImagesFromVideo.m on the alignment movies on both the .dat and the .avi files.

	Use createAlignment.m to make 3 alignment, follow the instructions to get the alignments in the correct order. Alignment files are saved in registration folder in 3dBrain, label them by date and the types of images aligned. The program will ask you to click matching points between the two images. 





STEP 1: WORM CENTERLINE DETECTION
(done locally for manual centerline initialization)

	initializeCLWorkspace.m - this will calculate background images and allow the user to manually initialize a few centerlines to help the detection algorithm.
 
				 - when this is done, move the alignment and the initialCLWorkspace.mat into the corresponding folder on tigress.
 
	Python submission code:

		submitWormAnalysisCenterline.py
	Matlab analysis code:

		clusterWormCenterline.m
	File Outputs:
	CLstartworkspace.mat, initialized points and background images for darkfield images
	CL_files folder, containing partial CL.mat files
	BehaviorAnalysis folder, containing the centerline.mat file with XY coordinates for each image.




	*NOTE: due to poor image quality of dark field images, it may be necessary to use some of the code developed by AL to manually adjust centerlines

STEP 2: STRAIGHTEN AND SEGMENTATION
	Python submission code:
		submitWormStraightening.py
	Matlab analysis code:
		clusterStraightenStart.m
		clusterWormStraightening.m
	File Outputs:	startWorkspace.mat, initial workspace used for during straightening for all volumes
			CLStraight* folder, folder containing all saved straightened tif files and results of segmentation.

STEP 3: NEURON REGISTRATION VECTOR ENCODING AND CLUSTERING
	Python submission code:
		submitWormAnalysisPipelineFull.py
	Matlab analysis code:
		clusterWormTracker.m
		clusterWormTrackCompiler.m
	File Outputs:	TrackMatrixFolder, containing all registrations of sample volumes with reference volumes.
			pointStats.mat, struccture containing all coordinates from all straightened volumes along with a trackIdx, the result of initial tracking of points. 

STEP 4: ERROR CORRECTION

	Python submission code:

		submitWormAnalysisPipelineFull.py

	Matlab analysis code:

		clusterBotChecker.m

		clusterBotCheckCompiler.m

	File Outputs:
 	botCheckFolder- folder containing all coordinate guesses for all times, one mat file for each neuron.
	pointStatsNew.mat- matfile containing the refined trackIdx after error correction. 

STEP 5: SIGNAL EXTRACTION

	Python submission code:

		submitWormAnalysisPipelineFull.py

	Matlab analysis code:

		fiducialCropper3.m

	File Output:
	heatData.mat, all signal results from extracting signal from the coordinates. 
	

#########################################################################
USEFUL GUIS FOR VISUALIZATION

ScanBinaryImageStack.m - Gui to view raw .dat file movies. This also works with .avi files
wormCLviewer.m - Gui to view darkfield worm images along with the centerline
VisualzeWorm3danalysis.m - Gui to view straightened worm data along with tracked coordinates
VisualizeTrackedData.m - Same as previous, but works on the unstraightened .dat file. 

#########################################################################
FOLDER DATA REQUIREMENTS

======Worm Videos======
sCMOS_Frames_U16_1024x1024.dat  -	binary image file from the Leifer Lab's WholeBrainImaging labview program. The images do not need to be 1024 by 1024.The image is produced by a dual view setup so that the left and right images represent the RFP and the GCaMP6s images in some order. 

LowMagBrain* folder containing all low magnification data including
cam0.avi	-	low magnification fluorescent images of the worm’s brain
cam1.avi	-	low magnification dark field images of the worm’s posture
CamData.txt	-	text file with relative timing for every frame

*****NOTE: FOR OLDER DATA
We used to use avi’s that were not time synced and used a YAML file containing the meta data. Using this data requires the code from the repo https://github.com/leiferlab/MindControlAccessUtils.git. 

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



