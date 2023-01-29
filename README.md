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
done, you can compile and test the algorithms by running this commands
```
mkdir build
cd build
cmake ..
make
```
you can build the project.
