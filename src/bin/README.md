
# DE424 MINI PROJECT: FIND THE WUMPUS

## Student Name: Amy McDermott
## Student Number: 26911264

### DIRECTORY STRUCTURE

The core files for this project are located as follows:

C++ Source Files:
- emdw/src/bin/wumpus1.cc (Dataset 1 code)
- emdw/src/bin/wumpus2.cc (Dataset 2 code)
- emdw/src/bin/wumpus3.cc (Dataset 3 code)

Visualisation:
- emdw/src/bin/datasetplot.py (Python visualisation tool)

Data Dependencies:
- Ground Truth: emdw/src/groundtruth/wumpus_trajectory[1-3].txt
- Observation Data:  emdw/src/dataset[1-3]/ (containing detection data)
- Model Output: emdw/build/wumpus_location[1-3].txt (Generated after run)

### COMPILATION AND EXECUTION (C++)

To compile and run the inference for any dataset (using Dataset 3 as an 
example), execute the following commands from the root 'emdw' directory:

    cd build
    cmake ../
    make -j7 wumpus3
    src/bin/wumpus3
    cd ..

Note: Replace 'wumpus3' with 'wumpus1' or 'wumpus2' for the 
respective datasets.

Upon successful execution, the inferred trajectory will be saved to:
emdw/build/wumpus_location[X].txt


### VISUALIZATION (PYTHON)

The Python script generates an animated GIF comparing the inferred 
trajectory (build folder) against the ground truth (groundtruth folder) 
and detections (dataset folders). In addition, the script prints a 
summary of quantitative metrics to the terminal so that model performance 
can be evaluated

Navigate to the bin directory and run the script with the dataset 
number as a command-line argument:

    cd src/bin
    python3 datasetplot.py [1/2/3]

Example for Dataset 3:
    python3 datasetplot.py 3

Output:
The resulting animation is saved as 'wumpus_tracking.gif' in the 
'emdw/src/bin' directory.
Accuracy, RMSE, and Mean Manhattan Distance (MMD) output to the terminal.


### REQUIREMENTS

- C++: emdw library dependencies (as per module setup).
- Python 3: requires matplotlib and numpy libraries.

--
### INSTALLATION AND SETUP FOR MARKING

To ensure the project compiles correctly within the standard DE424 
environment, please follow these placement instructions:

1. Copy the provided source files (wumpus1.cc, wumpus2.cc, wumpus3.cc) 
   into your existing 'emdw/src/bin/' directory.
2. Copy 'datasetplot.py' into the same 'emdw/src/bin/' directory.
3. Ensure the 'dataset1', 'dataset2', 'dataset3', and 'groundtruth' 
   folders are placed in 'emdw/src/' so that relative data paths 
   remain valid. The data provided from StemLearn must be in these 
   folders, named identically as they were provided.
4. Update 'emdw/src/bin/CMakeLists.txt' by adding the following lines 
   at the end:

   add_executable(wumpus1 wumpus1.cc)
   target_link_libraries(wumpus1 emdw)

   add_executable(wumpus2 wumpus2.cc)
   target_link_libraries(wumpus2 emdw)

   add_executable(wumpus3 wumpus3.cc)
   target_link_libraries(wumpus3 emdw)

Once these files are in place, proceed with the standard build 
instructions provided in Section 2.
