Python module for accumulating flow on a D8 flow grid. This class can be used to calculate drainage area and discharge,
and to accumulate any other tracer across a drainage network. 

Builds a network of nodes from a D8 flow grid. Uses a stack-based algorithm to traverse the network in topological order, 
modified from Braun & Willet (2013) DOI: 10.1016/j.geomorph.2012.10.008. This is faster than the recursive algorithm used in 
original Landlab implementation as we use an iterative build_stack algorithm (much faster). Most of the code is written 
in Cython for speed. The approach is linear w.r.t. the number of nodes in the network. Class is designed to be used with 
geospatial rasters, but can also be used with a numpy array of D8 flow directions. 

## Installation 

To install run: 
```
python setup.py build_ext --inplace
```
This compiles the Cython scripts and allows them to be imported into python scripts. Note this does not install a package so 
any scripts using this module must be in the same directory as the compiled files.

## Example of use 

An example of use is given in the script `example.py` which analyses the flow accumulation and drainage area of the 
domain covered by `d8_eg.nc`. 