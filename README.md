# isostasy
Stand-alone isostasy model

This compact module is built to run simple bedrock isostasy models.
So far, the implementations include:

1. Local lithosphere, relaxing asthenosphere (LLRA)
2. Elastic lithosphere, relaxing asthenosphere (ELRA)
3. Elementary GIA model (EGIA) - work in progress!

## Initial configuration

1. Define a configuration file for your system and compiler within the folder `config`. Use a previously existing config file as a template. 

2. Generate a Makefile using the `config.py` script:

```
python config.py config/pagos_gfortran
```

where the path to the config file should match your desired choice. At this point, you should have a valid Makefile in the main folder that can be used to compile the program. 

3. Compile the test program

```
make test_isostasy
```

4. Run the test program

```
./lib/bin/test_isostasy.x 
```

The output directory is hard-coded to `output/test-isostasy` for now. When the program runs successfully, the model output will be saved in the file `bedtest.nc` in the output directory. To play with different parameters, modify them in the default parameter file `par/test_isostasy.nml`. 