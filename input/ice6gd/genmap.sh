# usage example: ./genmap.sh (arguments provided inside)

domain=Antarctica

grid_name_src=ice6gd     # Source grid (input file: grid_ice6gd.txt)
grid_name_tgt=ANT-32KM   # Target grid (input file: grid_ANT-32KM.txt); possible values: ANT-4KM, ANT-8KM, ANT-20KM, ANT-32KM

#nc_src=$home/models/EURICE/ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 

nc_src=jgrb52450-sup-0005-data_s5.nc  # Input file contaning data to be gridded ; output will be e.g. scrip-con_ice6gd_ANT-32KM.nc

cdo gencon,grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc


# To perform remapping using the weights file
# cdo remap,grid_${grid_name_tgt}.txt,scrip-con_${grid_name_src}_${grid_name_tgt}.nc ${infile} ${outfile}

