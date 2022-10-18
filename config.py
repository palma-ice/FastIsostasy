'''
Script to generate Makefile with the correct compiler
configuration for a given machine.
'''

from subprocess import *
import argparse
import sys 

# Manage command-line arguments to load compiler option
parser = argparse.ArgumentParser()
parser.add_argument(
    'config', metavar='CONFIG', type=str,
     help='Name of config file to use (e.g., config/pagos_gfortran)')
args = parser.parse_args()


# Determine locations
target_dir  = "./"
config_dir  = "config/"
config_path = args.config 

# Load template Makefile and insert compiler configuration
makefile0    = open(config_dir+"Makefile").read()
compile_info = open(config_path).read()
makefile1 = makefile0.replace("<COMPILER_CONFIGURATION>",compile_info)

# Write the new Makefile to the main directory
open(target_dir+"Makefile","w").write(makefile1)

print( "".join(["\n\nMakefile configuration complete for configuration file: ",config_path,"\n"]) )

