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



