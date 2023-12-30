# Perform remapping using cdo by defining the source and target grid_name values,
# as well as the source file infile and the target file outfile. This script
# assumes that the scrip map weights have already been generated, e.g. using genmap.sh.

# usage example: ./remap.sh ice6gd ANT-32KM jgrb52450-sup-0005-data_s5.nc H_ice6gd_32KM.nc

grid_name_src=$1
grid_name_tgt=$2
infile=$3
outfile=$4

# To perform remapping using the weights file
cdo remap,grid_${grid_name_tgt}.txt,scrip-con_${grid_name_src}_${grid_name_tgt}.nc ${infile} ${outfile}
