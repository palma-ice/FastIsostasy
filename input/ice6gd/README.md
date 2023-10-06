
To remap a data file, 

1. Define the target and source grids using the format required by `cdo`. If a grid is a simple lat-lon grid, or the projection is well, defined, this can be done via `cdo griddes file.nc > grid_file.txt`. 

2. Run the `genmap.sh` script, defining the options at the beginning of the script. Essentially this script just calls the cdo command `cdo gencon` to produce conservative mapping weights for going from the source to target grid.

3. Run the `remap.sh` script, to actually remap a data file defined on the grids of interest. Make sure to define the input arguments corretly. Again, this script essentially just calls `cdo remap` using the weights defined earlier.

That's it!
