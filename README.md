# Array Designer #
### Created by Dr. Rob Cooper and Dr. Sabrina Brigadoi ###

#########################################################
If you use any of these tools please cite:

Sabrina Brigadoi, Domenico Salvagnin, Matteo Fischetti, Robert J. Cooper, "Array Designer: automated optimized array design for functional near-infrared spectroscopy," Neurophotonics 5(3) 035010 (13 September 2018) https://doi.org/10.1117/1.NPh.5.3.035010
#########################################################

### Introduction ###:
This repository is a set of GUIs and scripts that create an easy-to-use environment to run the automated fNIRS array design methodology described by Brigadoi et al. in the paper listed above.

The repository consists of two matlab apps (ArrayDesignerAPP and ROIBuilderAPP) and a range of directories containing associated scripts and data.

The basic premise of ArrayDesigner (AD) is that a user can select the areas of a cortical mesh model they wish to target with their fNIRS array, then input parameters associated with their hardware (e.g. number of source and detector locations) and AD will automatically generate optimised positions for those sources and detectors on the scalp (in the 10-5 or 10-2.5 system) to maximise the sensitivity (or sensitivity and coverage) to/of that target cortical area.

### Installation ###:

1) To run Array Designer you need Matlab - ideally 2019b or later. Depending on the version you download, you will obtain a main directory ("ArrayDesigner") containing several subdirectories and .mlapp files. This entire directory simply needs to be added to you matlab path.

2) If you downloaded a separate PMDFs.data file, this needs to be placed within the directory of AD's default head model at: ArrayDesigner/HeadModels/MNI152_Symmetric/PMDFs/PMDFs.data

3) If you are using a Mac, you will need to bypass the security settings that prevent the running programs from unknown sources. In the Finder on your Mac, locate the ArrayDesigner/Build/nirsmain file. Control-click the app icon, then choose Open from the shortcut menu. This won't run properly, but the app will be saved as an exception to your security settings, and ArrayDesigner will then be able to ope it when it needs it.

4) If you want to use a different head model to the default (which we encourage you to try!) you will need TOAST++ downloaded and in your matlab path.(http://web4.cs.ucl.ac.uk/research/vis/toast/).

### ArrayDesignerAPP ###:
To run ArrayDesignerApp.mlapp just double click on it. Note if you double click on it via matlab, the App Designer interface will open - you will need to click 'Run' to run it. The main app GUI is arranged in hierarchical/logical order from top to bottom. In order of operation, users should
1) Select a head model in which to operate
2) Select an ROI to target
3) Select a solution space (10-5 or 10-2.5 (default))
4) Enter the number of sources and detectors
5) Set the coverage weight parameter
6) Set the source-detector/optode separation parameters
7) RUN ARRAY DESIGNER

Details:
On (1) - The default head model is based on the the MNI152 symmetric atlas. This and any other head models are stored in the AD repo /HeadModels. Each model has a directory, within which is a .mshs file. The .mshs format (which is just a matlab .mat) is described in the DOT-HUB toolbox. It contains a volume and GM surface mesh and various other variables. AD can work with any head model in the .mshs format, so if you want to move beyond the atlases provided, all you need to do is create a .mshs file for your head meshes, save it in the appropriate folder (e.g. AD/HeadModels/myMesh/myMesh.mshs) and load it. Note that we haven't had much time to test this yet, but there is no reason it should not work. Note that GM surfaces should be relatively sparse (circa 10,000 nodes) to maintain speed, and parcellations of the GM surface (included in the variable structure gmSurfaceMesh.parcellation) are useful for building ROIs.

Within the HeadModels folder for the MNI_Symmetric model (which is the default for AD) are various other things used by AD including /PMDFs/PMDFs.data. This is a large (~1Gb) binary containing all possible PMDFs for the 10-2.5 solution space for the MNI_Symmetric model. This allows users to run AD using the default head model without having to calculating PMDFs themselves. If you select the MNI_Asymmetric model (or any other), you will have to calculate the PMDFs locally. (We didn't want to include 2x1Gb binary files in the release, but the asymmetric model is the one used in the paper, so we thought it important to include it - the ROIs from the paper are also in that folder). AD does let you calculate PMDFs, and the process should be straightforward (if a little time consuming), but you will need to have TOAST++ installed.

On (2) - The ROI is a .txt file that represents the locations on the selected head model's GM surface mesh that you wish to target with your array. The ROIs should be saved in the relevant head model folder in a subdirectory named 'ROIs' (e.g. HeadModels/MNI_Symmetric/ROIs/myROI.txt). To make an ROI, go to the Tools menu, and open ROIbuilder.

On (5) - The objective function that Array Designers uses is essentially = [total sensitivity + Cw * coverage] where Cw = coverage weight. If Cw = zero, Array Designer will return the position of sources and detectors on the scalp that maximise sensitivity to the target cortical area, but unless your ROI is very focal, this is rarely what you actually want from an fNIRS array. Increase the coverage weight to force the array to spread out over your ROI.

On (6) All these parameters are as defined in the paper, but they should hopefully be self explanatory.

Note that the inputs and results of a run of AD are saved as .AD files in /Results whenever you run Array Designer. You can load prior .AD files to save you having to re-enter the same information again etc.


### ROIBuilderAPP ###:
This app allows you to create an ROI - which specifies the area of the cortex you want your array to be optimised to cover.

Once a head model is selected, a cortical surface is displayed. If that head model includes parcellations (as the MNI_Symmetric does), the parcel list view can be used to select parcels to add to or remove from your ROI. You can also use the Brain Painter controls section to paint areas of your ROI. Simply click the data pointer onto a location on the cortex figure, set your brush diameter and hit 'paint'.

Once complete, you can save your ROI .txt file. Remember this needs to be saved in the relevant head model folder in a subdirectory named 'ROIs', (e.g. HeadModels/MNI_Symmetric/ROIs/myROI.txt).

Happy Designing!


