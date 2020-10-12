# Aneurysm Thickness Measurement

# General Info
This program was designed to take an input of 2 .off files, convert that data into vtkPolyData, and compute thickness measurements between the two geometric shapes. This implementation is specifically designed for measuring the thrombus thickness of abdominal aortic aneurysms.

# Compilation
The CMake file included appropriately compiles all necessary code with the assumption that the system has Visualization Toolkit installed on their system. You can compile the code by running

```
make
```
# Running the Program
To start the program, simply run the command listed below. Once started, the program will prompt a user for a file of the outer aorta wall. Then the program will prompt a user for a file representing the lumen wall of the aorta. These must both be files of type .off.

```
./GeometryModule
```

# Testing 
Excess code which was not utilized exemplifies the various methods attempted in creating precise, reproducible results. The final functionality determined to be the most accurate was an inner-to-outer measurement of the nearest neighbouring node on the outer mesh for every point on the inner mesh. 

The optimal point density was also researched in testing. The section to accumulate results with varying point density can be uncommented to determine the optimal point density for specific meshes. A point density of 0.328 points/mm^2 was determined to balance accuracy and computation time on average.

# Saved Results
The results of the nearest neighbour algorithm from outer-to-inner are saved in a .vtp file titled "wall.vtp". The results of the nearest neighbour algorithm from inner-to-outer are saved in a .vtp file titled "lumen.vtp". 

