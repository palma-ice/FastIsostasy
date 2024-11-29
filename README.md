# FastIsostasy

FastIsostasy is a model that regionally computes the glacial isostatic adjustment (GIA), as described in [Swierczek-Jereczek et al. (2024)](https://egusphere.copernicus.org/preprints/2023/egusphere-2023-2869/#discussion). It approximates the gravitational response, accounts for the resulting heterogeneity of the sea-surface height and computes the evolution of masks (continent, ocean, floating ice, grounded ice), as well as the load anomalies applied upon the solid Earth. The resulting bedrock deformation can be computed with different models:

1. Local lithosphere, relaxed asthenosphere (LLRA).
2. Elastic lithosphere, relaxed asthenosphere (ELRA; LeMeur and Huybrechts, 1996).
3. Elastic lithosphere, relaxed mAntle (ELVA; Cathles, 1975; Lingle and Clark, 1985; Bueler et al., 2007).
4. Laterally-Variable ELVA (LV-ELVA; Swierczek-Jereczek et al., 2024).

Remark: LV-ELVA is a generalisation of ELVA and both therefore rely on the same backend. To use ELVA, simply pass laterally-constant solid-Earth parameters to the solver.

If you are interested in coupling FastIsostasy to your ice-sheet model, the remainder of the README.md is for you! It provides you with the initial steps to get started with FastIsostasy (stand-alone) and the few lines of fortran code needed for a correct interaction with your ice sheet.

If compatible with your research ecosystem, you should preferably use the [Julia implementation](https://github.com/JanJereczek/FastIsostasy.jl), which offers extend capabilities (GPU acceleration, higher-order time integration schemes, adaptive time stepping, time-variable ocean surface). 

## Getting started (stand-alone)

1. Define a configuration file for your system and compiler within the folder `config`.
    Use a previously existing config file as a template. 

2. Generate a Makefile using the `config.py` script:

    ```
    python config.py config/pik_hpc2024_ifx
    ```

    where the path to the config file should match your desired choice. At this point, you should have a valid Makefile in the main folder that can be used to compile the program. 

3. Compile the test program

    ```
    make test_isostasy
    ```

    You can choose which test should be performed by setting `experiment` (defined in `test_isostasy.f90`) accordingly.

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

Coupling FastIsostasy requires three important steps:
- initialize the module on the desired domain,
- initialize the state of FastIsostasy,
- update the state of FastIsostasy.

The fortran code corresponding to these steps is (`! ...` refers to segments of your own code):

```
! ...

type(isos_class)        :: isos1
character(len=512)      :: path_par
character(len=256)      :: file_isos_restart

call isos_init(isos1, path_par, "isos", nx, ny, dx, dy)

! ...

call isos_init_state(isos1, z_bed, H_ice, time)

! ...

do time = 1, nt
    ! ...
    call isos_update(isos1, H_ice, time, dwdt_corr=dzbdt_corr)
    z_bed = isos1%out%z_bed
    z_sl  = isos1%out%z_ss
    ! ...
end do

call isos_restart_write(isos1, file_isos_restart, time)

! ...
```