# CLG

A C++ implementation of the Conditional Linear Gaussian model and its inference
methods such as they are presented in "Graphical Models for Associations between
Variables, some of which are Qualitative and some Quantitative by S. Lauritzen and
N. Wermuth (1989)" and "Propagation of Probabilities, Means, and Variances in Mixed
Graphical Association Models by S. Lauritzen (1992).". 

WARNING: At the moment, the marginalization of discrete variable is not implemented
and the method can only be used for continuous only random variables.

This implementation uses the C++ libraries aGrUM (for graphical models) and Eigen
(for matrix operations) and you need to install them before using it. Once it is
done, you can compile and test the algorithms by first modifying the 
lines:
```
set(AGRUM_INSTALLATION_DIRECTORY "/home/lasserre/.venv/otagrum")
set(aGrUM_DIR "${AGRUM_INSTALLATION_DIRECTORY}/home/lasserre/.venv/otagrum/lib/cmake/aGrUM/")
```
in the CMakeList.txt file to specify the location of aGrUM on your computer.
Then, by running this commands in the cpp/src directory:
```
mkdir build
cd build
cmake ..
make
```
you can build the project.
