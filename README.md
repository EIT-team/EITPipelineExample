# Making EIT images
How to go from Geometry -> Mesh -> Forward -> Inverse -> EIT image

To make an image you need to install and make visible to Matlab:

0. [iso2mesh](https://github.com/fangq/iso2mesh) - add folder to matlab path. This contains useful tools for handling meshes in Matlab
0. [Mesher](https://github.com/EIT-team/Mesher) - add MATLAB folder to path. If you are using WSL to run the linux software on Windows, its better to have this in a windows folder that in ubuntu looks like `/mnt/c/Users/User/Mesher`. This creates the FEM models from segmentations
0. [Supersolver](https://github.com/EIT-team/SuperSolver) - add src directory. This runs the forward model and creates the jacobian
0. [Reconstruction](https://github.com/EIT-team/Reconstruction) - add src/matlab. This inverts the jacobian and creates images
0. [paraview](https://www.paraview.org/download/) - views the output vtk files

Each directory has examples, some of which is repeated here, so if something is not working then check within the specific directory for more info. 