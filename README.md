# Nucleolus_simulation
These codes are all for the simulation of sub-nucleolus simulation

Mento: Searching for the optimal distribution of processing factors in one sub-nucleolus unit using Mento Carlo strategy

SimulationClusterNumber: Simulation the cluster number in the max-cross section from different cluster number of one sub-nucleolus unit

ImageAlignment: Register two images to get an average image (This code is adapted from the https://ww2.mathworks.cn/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation)

## Mento
inner3_calo_meto_v3.m: The main matlab code for simulation. In a summary, all the processing factors are firstly randomly distributied in a 50*50*50 3D space. Then calculate the binding ability of each processing factor and give a disturbition to the 500 processing factors which have the lowest binding ability to pre-rRNAs, which are moved to the positions with a higher binding ability. After giving 10000 disturbition to this system and get the finnally optimal distribution of processing factors in a sub-nucleolus region 

moveSteps_max_v1.m: This function find the 500 processing factors with the lower binding ability and randomly decide the moving direction and steps of them to new position of the higher binding ability 

CalEffect_all_new_v4.m: Code for calculation of binding ability for each processing factor. The formula for calculation of binding ability of each processing factor is f(x,y,z)=p(x,y,z)/sum(p(delta(l)*delta(l)). (x,y,z) is the coordinate of processing factor, p(x,y,z) is the processing factor concentration at coordinate (x,y,z), and l is the path from the FC/DFC border to the outside surface of the DFC.

visualization_Mento_Carlo.m: Code for visualization of mento_carlo simulation process.

random_v4_s3_4_500_500_res.mat: the matrix record each disturbution, the results from inner3_calo_meto_v3.m and the input for the visualization_Mento_Carlo.m

change_color.avi: The video for visualization simulation process and it is the result from the visualization_Mento_Carlo.m

### Usage
Directly run inner3_calo_meto_v3.m in Matlab, then get the matrix random_v4_s3_4_500_500_res.mat, finnaly run the code visulization_Mento_Carlo.m and get the visualization result for Mento Carlo simulation

## SimulationClusterNumber
cutSphere.m: this code is the main for the cluster number simulation. A given number (12-30) of clusters (solid sphere) are individually arraged around the original point in the shell of DFC sphere with the radius of 247.5 nm, which was measured under SIM images. The position of clusters was randomly chosen for 100 times, images of the max-cross section were collected from 200 different angles randomly selected. This code automately count the number of cluster in the max-cross section. The output files are 100 text file record the number of clusters in max-cross section in each randomly chosen clusters generating.

cutSphere_visulization_version.m: This code is for visulization of the image of max-cross section from 200 different angles of one randomly chosen number of cluster. The output files are 200 tiff files of the max-cross section images.

rotateEllipsoid.m: One function for rotating the clusters in the center of origin point.

### Usage
Change the [number_shape] in the cutSphere.m from 12 to 30, then run this code in matlab and there will be 100 text files recording the the number distribution of max-cross image.
Change the [number_shape] in the cutSphere.m from 12 to 30, and will get 200 tiff files for the images of max-cross section in a given number of clusters (12-30).

## ImageAlignmet
This code is adapted from https://ww2.mathworks.cn/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) for image registration and alignment

FBL_alignment.m: Main code for image alignmet. Two input images, one is as a reference image and another is the aligning image. for the aligning image, it is rotated from 0 to 359 degree by 1 degree per rotation and the rotated images are stored after each rotation. After applying fast Fourier transform to the reference image and rotated images, the cross-correlation peak was individually calculated between the reference image and each rotated image. A rotated image, which has the best cross-correlation peak with the reference image, is then refined in the x-y scale according to the reference image to be the registered image. And 


