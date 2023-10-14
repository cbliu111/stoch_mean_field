### Stochastic simulation of chemical reaction network with mean-field randomization.

The repository contains the codes essential for generating all the results presented in the paper "Error-controlled coarse-graining dynamics with mean-field randomization". The implementation here provides a C++ implementation of all the algorithms included in the paper, including RNRM, R-leap and $\omega$-leap. All these algorithms are based on the mean-field randomization method, and provide coarse-grained sampling of trajectories of chemical reaction networks. The application of these algorithms are straightforward by following the discussions in the paper. 

### Building

The building process for these algorithms involves using cmake to generate the necessary makefiles for building the C++ files. To build the project, follow these steps:
1. Clone the repository to your local machine.
2. Install cmake if it is not already installed on your system.
3. Create a build directory within the repository directory.
4. Navigate to the build directory and run "cmake .." to generate the makefiles.
5. Run "make" to build the project. 
Once the project is built, you can run the executable to simulate the behavior of the reaction process in the model. The parameters for the simulation can be adjusted in the source code. 

Note the $\omega$-leap in the paper is originally named as the "adaptleap" algorithm. 

