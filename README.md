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

### ArrayDesignerAPP ###:
The main app is arranged in hierarchical/logical order from top to bottom. In order of operation, users should
1) Select a head model in which to operate
2) Select an ROI to target
3) Select a solution space (10-5 or 10-2.5 (default))
4) Enter the number of sources and detectors
5) Set the coverage weight parameter
6) Set the source-detector/optode separation parameters
7) RUN ARRAY DESIGNER

Details:
On (1) - The head models are stored in the AD repo /HeadModels. Each model has a directory. Two adult atlas head models are provided with AD. MNI_Symmetric and MNI_Asymmetric. The model itself is a .mshs format (which is just a matlab .mat) as described in the DOT-HUB toolbox. It contains a volume and GM surface mesh and various other variables. AD can work with any head model in this format, so if you want to move beyond the atlases provided, all you need to do is create a .mshs file for your head meshes, save it in AD/HeadModels/ and load it. Note that we haven't had much time to test this yet, but there is no reason it should not work. Note that GM surfaces should be relatively sparse (circa 10,000 nodes) to maintain speed. 

Within the HeadModels folder for the MNI_Symmetric model (which should be considered the default for AD) are various other things used by AD including /PMDFs/PMDFs.data. This is a large (~1Gb) binary containing all possible PMDFs for the 10-2.5 solution space for the MNI_Symmetric model. This allows you to use this model without calculating PMDFs yourself. If you select the MNI_Asymmetric model, you will have to calculate the PMDFs locally. (We didn't want to include 2x1Gb binary files in the release, but the asylum model is the one used in the paper, so we thought it important to include it - the ROIs from the paper are also in that folder). To calculate PMDFs, you will need to have TOAST++ installed.

On (2) - The ROI is a .txt file that represents the locations on the selected head model's GM surface mesh that you wish to target with your array. The ROIs should be saved in the relevant head model folder in a subdirectory named 'ROIs' (e.g. HeadModels/MNI_Symmetric/ROIs/myROI.txt). To make an ROI, go to the Tools menu, and open ROIbuilder.

On (5) - The objective function that Array Designers uses is essentially = [total sensitivity + Cw * coverage] where Cw = coverage weight. If Cw = zero, Array Designer will return the position of sources and detectors on the scalp that maximise sensitivity to the target cortical area, but unless your ROI is very focal, this is rarely what you actually want from an fNIRS array. Increase the coverage weight to force the array to spread out over your ROI.

On (6) All these parameters are as defined in the paper, but they should hopefully be self explanatory.

Note that the inputs and results of a run of AD are saved as .AD files in /Results whenever you run Array Designer. You can load prior .AD files to save you having to re-enter the same information again etc.



### ROIBuilderAPP ###:
This app allows you to create an ROI - which specifies the area of the cortex you want your array to be optimised to cover.

Once a head model is selected, a cortical surface is displayed. If that head model includes parcellations (as the MNI_Symmetric does), the parcel list view can be used to select parcels to add to or remove from your ROI. You can also use the Brain Painter controls section to paint areas of your ROI. Simply click the data pointer onto a location on the cortex figure, set your brush diameter and hit 'paint'.

Once complete, you can save your ROI .txt file. Remember this needs to be saved in the relevant head model folder in a subdirectory named 'ROIs', (e.g. HeadModels/MNI_Symmetric/ROIs/myROI.txt).

Happy Designing!


