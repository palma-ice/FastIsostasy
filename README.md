# isostasy
Stand-alone model of regional glacial isostatic adjustment (GIA).

This compact module is built to run regional GIA models, which compute the deformational
response to changes in the surface load and, optionally, the gravitational response,
the heterogeneous sea-surface height and the sea-level contribution, as defined in
Goelzer et al. (2020) and extended in Swierczek-Jereczek et al. (2024).

So far, the implementations of the deformational response include:

1. Local lithosphere, relaxing asthenosphere (LLRA).
2. Elastic lithosphere, relaxing asthenosphere (ELRA). By providing a heterogeneous field
 of the relaxation time, the laterally-variable formulation of Coulon et al. (2021) is used.
3. Laterally-Variable Elastic Lithosphere/Viscous Asthenosphere (LV-ELVA), as described
   in Swierczek-Jereczek et al. (2024). Note that this is a generalisation of ELVA, as
   described in Bueler et al. (2007).

## Initial configuration and running the program

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


That's it!