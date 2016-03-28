 
Worm analysis protocol:


This repository hold the code used for the analyzing movies from the Leifer Lab's Whole brain imaging set. 




Folder requirements:
sCMOS_Frames_U16_1024x1024.dat  -	binary image file from the Leifer Lab's WholeBrainImaging labview program. The image is not actuallylly 1024x1024, but 1200x600 unsigned 16 bit integers.


LabJackData.txt - tab separated data file from lab jack from Leifer Lab's WholeBrainImaging labview program. Has stage, function generator, and Piezo data, recorded at 1kHz. 


CameraFrameData.txt - tab separated data file from lab jack from Leifer Lab's WholeBrainImaging labview program. Has image data, and is output at the rate of image capture (200Hz)



*fluor*.avi -	 video of low magnification fluorescent video from the Leifer lab's Colbert code 

*fluor*.yaml -	 metadata for low magnification fluorescent video from the Leifer lab's Colbert code, includes timing data

*behav*.avi - 		video of low magnification behavior video from the Leifer lab's Colbert code, will be used for centerline extraction 

*behav*.yaml -	 metadata for low magnification behavior video from the Leifer lab's Colbert code, includes timing data

alignments.mat -matlab structure with fields containing the affine transforms for image alignment. Alignments created using the createAlignment.m program: The file has the following fields
lowResFluor2BF	- affine transform between low magnification fluorescent and behavior images
S2AHiRes	- affine transformation between the RFP channel and the Gcamp6s channel of the hi magnification image along the the cropping rectangles for each of them.
Hi2LowResFluor	- affine transformation between the RFP hi magnification image (before cropping) and the hi magnification fluorescent images






#########################################################################

Main scripts in analysis pipeline:

STEP 1: CENTERLINE DETECTION 
tripleImageCenterline_par.m

STEP 2: STRAIGHTEN AND SEGMENTATION 
StraightenWormFolder20150910.m
STEP 3: TRACKING ON CLUSTER USING PYTHON SCRIPTS submitWormAnalysisPipeline.py
STEP 4: COMPILE CLUSTER RESULTS AND FILL HOLES compileBotChecker.m
STEP 5: CROP REGIONS AROUND TRACKED POINTS 
fiducialCropper2.m




#########################################################################
Useful GUI programs
ScanBinaryImageStack.m - Gui to view raw .dat file movies. This alsow works with .avi files
wormCLviewer.m - Gui to view darkfield worm images along with the centerline
VisualzeWorm3danalysis.m - Gui to view straightened worm data along with tracked coordinates