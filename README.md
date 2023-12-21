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
This compiles the Cython scripts.

## Example of use 

```
from d8_accumulator import D8Accumulator, write_geotiff
accumulator = D8Accumulator("d8.tif")
# Create an array of cell areas
cell_area = np.ones(len(accumulator.receivers)) * 100 # 100 m^2 cell area
# Calculate drainage area in m^2
drainage_area = accumulator.accumulate(weights=cell_area)
# Calculate number of upstream nodes
number_nodes = accumulator.accumulate()
# Write the results to a geotiff
write_geotiff("drainage_area.tif", drainage_area, accumulator.ds)
```
