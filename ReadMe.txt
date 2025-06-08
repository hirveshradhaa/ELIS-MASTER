
# ELIS-MASTER

A C++ implementation of a time-series classification framework based on shapelet discovery and classification, designed to work with datasets from the [UCR Time Series Classification Archive](http://www.cs.ucr.edu/~eamonn/time_series_data/).

## Project Overview

ELIS-MASTER is a C++ project that implements a time-series classification framework using shapelets. Shapelets are small, discriminative subsequences within time series data that are highly effective in distinguishing between classes.

 What It Does

- Shapelet Discovery: Identifies the most relevant subsequences (shapelets) from the training dataset that best separate different classes.
- Shapelet Adjustment & Classifier Training: Fine-tunes shapelets to improve classification performance and train the model.
- Prediction: Uses the trained shapelets to classify unseen time series in the test set and output results.

Why It’s Important

- Supports **interpretable time-series classification**, crucial in fields like finance, healthcare, and IoT.
- Efficient and scalable implementation using C++ and standard build tools like CMake.
- Utilizes datasets from the [UCR Time Series Benchmark](http://www.cs.ucr.edu/~eamonn/time_series_data/) for standardized benchmarking.

## Key Features

- Discover shapelet candidates from training time series data
- Adjust shapelets and train the classifier
- Predict and classify test time series
- Modular and extensible C++ codebase



After downloading the file, all the required files in the Final Project Folder.

## Building the Project:

To construct the project, follow the steps outlined below within the src directory:

1. Create the run directory by executing the following command:
  
   mkdir run :—> make run directory
  
2. Inside the run directory, create a data subdirectory and copy-paste the dataset from the [UCR Benchmark](http://www.cs.ucr.edu/~eamonn/time_series_data/) as specified in the paper:
  
   mkdir run/data :—> here we have to copy paste the dataset from dataset in [UCR Benchmark](http://www.cs.ucr.edu/~eamonn/time_series_data/) as per mentioned in the paper

3. Execute the following command to run CMake in the current directory. This command is commonly used in C++ to execute CMake, a cross-platform build system:

   cmake . :—> The cmake . command is used in C++ to execute CMake, which is a cross-platform build system. When executed in a directory containing a CMakeLists.txt file, the cmake . command generates the build system for the project in that directory.

4. Use the make command in C++ to compile and build C++ programs based on the instructions provided in the makefile. A makefile is a script specifying how the program should be built, including the source files, dependencies, and compilation instructions:

   make :—> The make command in C++ is used to compile and build C++ programs based on the instructions provided in the makefile. A makefile is a script that specifies how the program should be built, including the source files, dependencies, and compilation instructions.








To store the output in a designated folder, follow these commands. As an example, if generating shapelets for the StarLightCurves dataset:

cd run
mkdir StarLightCurves
cd StarLightCurves

## How to Use:

1. To identify shapelet candidates from training time series, execute the following command:
  
   ../discover 20 ../data/StarLightCurves/StarLightCurves_TRAIN 30 > dout
 

2. To adjust shapelets and train the classifier, execute the following command:
  
   ../adjust train 20 ../data/StarLightCurves/StarLightCurves_TRAIN 30 -25 0.01 0 0.01 > tout
  

3. Finally, to predict for test time series, use the following command:

   ../adjust test 20 ../data/StarLightCurves/StarLightCurves_TEST 30 -25 > rout


The approach generates three text files:

- init.txt: Discovered shapelet candidates for each class
- learned.txt: Shapelets for each class
- result.txt: Testing results
