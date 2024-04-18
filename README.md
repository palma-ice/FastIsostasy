# FastIsostasy

FastIsostasy is a model that regionally computes the glacial isostatic adjustment (GIA), as described in [Swierczek-Jereczek et al. (2024)](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-2869/#discussion). It is designed to be easily coupled to an ice-sheet model and offers a [Julia implementation](https://github.com/JanJereczek/FastIsostasy.jl), in addition to the Fortran one hosted in the present repository. FastIsostasy approximates the gravitational response and accounts for the resulting heterogeneity of the sea-surface height when computing the evolution of the masks (continent, ocean, floating ice, grounded ice), as well as the load applied upon the solid Earth. The resulting deformation of the bedrock can be computed with different models:

1. Local lithosphere, relaxed asthenosphere (LLRA).
2. Elastic lithosphere, relaxed asthenosphere (ELRA; LeMeur and Huybrechts, 1996).
3. Elastic lithosphere, relaxed mAntle (ELVA; Cathles, 1975; Lingle and Clark, 1985; Bueler et al., 2007).
4. Laterally-Variable ELVA (LV-ELVA; Swierczek-Jereczek et al., 2024).

Remark: LV-ELVA is a generalisation of ELVA and both therefore rely on the same backend. To use ELVA, simply pass laterally-constant solid-Earth parameters to the solver.

## Initial configuration and running the program

1. Define a configuration file for your system and compiler within the folder `config`.
    Use a previously existing config file as a template. 

2. Generate a Makefile using the `config.py` script:

    ```
    python config.py config/pagos_gfortran
    ```

    where the path to the config file should match your desired choice. At this point, you should have a valid Makefile in the main folder that can be used to compile the program. 

3. Compile the test program

    ```
    make test_isostasy
    ```

4. (Optional) If you want to run more complex simulations (e.g. Test 4), you will need additional data (ice loading history, parameter fields, etc.), which can be downloaded via:

    ```
    git clone https://github.com/JanJereczek/isostasy_data.git
    git clone git@github.com:JanJereczek/isostasy_data.git      # if you are using ssh
    ```

    If you already downloaded `isostasy_data`, make a symbolic link to it:

    ```
    ln -s path_to/isostasy_data
    ```

5. Run the test program

    ```
    ./libisostasy/bin/test_isostasy.x 
    ```

The output directory is hard-coded to `output/` for now. When the program runs successfully, the model output will be saved in the file `output/bedtest_x.nc`, with `x` the experiment number associated with the parameters `par/test_isostasy_testx.nml`. The experiment number `x` can be changed within `tests/test_isostasy.f90`.

That's it!

## Coupling FastIsostasy to your favorite ice-sheet model

The interface of FastIsostasy is apparent in `tests/test_isostasy.f90` and we refer to the structure used there for any coupling to an ice-sheet model. For Julia users, please consult [this page](https://janjereczek.github.io/FastIsostasy.jl/dev/examples/tutorial/#Make-your-own-time-loop).