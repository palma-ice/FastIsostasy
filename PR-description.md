# Refactoring for modularity, performance and sea-level treatment

This PR contains a substantial refactoring of the isostasy module. It does not comply
with the guideline of atomic PRs (one PR = one task) because of the time constraints.
However, all tests are passed by the refactor and a detailed description follows.

## Major changes

Modularity:
- `isostasy.f90` has been splitted into `convolutions.f90`, `finite_differences.f90`,
 `green_functions.f90`, `isos_utils.f90`, `kelvin_function.f90`, `sealevel.f90` and
 `lv_xlra.f90`. The latter includes the functions solving the viscous displacement
 of ELRA and LLRA (gathered under XLRA), optionally with laterally variable relaxation
 time scale `tau`.
- `solver_elva.f90` has been merged with `lv_elva.f90` since it avoids duplication.
 The code is now sufficiently performant for LV-ELVA to eb just as quick as ELVA.
- Along these lines, there are now only 4 methods that can be used to solve the viscous
 displacement (0: constant displacement, 1: LLRA, 2: LV-ELRA, 3: LV-ELVA).
- `isostasy_defs.f90` defines a new derived type `isos_domain_class`, that stores all
 the information related to the computation domain (including the pseudo-differential
 operator in Fourier space as well as the FFT plans).
- `isos_class` now has `ref` as a field, which contains the values of `now` at
 initialisation.
- Many functions had duplicates, which has been corrected. Example: function computing
 `kappa`, the pseudo-differential operator in Fourier space.
- The code that was used several times has now been written as functions. Example:
 computing the rigidity of the elastic lithosphere based on other parameters.


Performance and accuracy:
- The normalisation of the FFT is only applied on the backtransform, as it should be!
- Convolutions are computed based on FFTs.
- FFT plans are precomputed.
- FFTs of Green's functions are precomputed (elastic, viscous and ssh).
- The convolution of the viscous Green function with the load gives the ELRA equilibrium
 displacement. This solves the trade-off between accuracy and performance that previously
 existed because of the time-domain convolution. In particular, we can just use
 `radius_fac = 6.0` and `filter_scaling = 1.0` (more accurate than the pre-existing scaled
 version) without slowing the code.
- Allocations inside the time loop are reduced. As Alex suggested, this only improves
 the computation time marginally.

Regional sea-level model (ReSeLeM):
- Column anomalies of ice, seawater, lithosphere and mantle are now used as forcing.
- "Geoid" (somewhat obsolete word) replaced by sea-surface height and its perturbation.
- Introduction of coherent masks (ocean, continent, grounded and active region).
- ReSeLeM + LV-ELRA gives elementary GIA model of Coulon et al. (2021)

General physics:
- Allow rectangular domain (without the need of a square extension). The function that
 treats the problem by extending the fields on a square domain is however still here
 for legacy purposes. Sometimes, an extension of the domain is needed but it should
 have a similar x-y aspect ration than the original domain. This is a task for the future!
- Distortion matrix included in computations. However, the matrix itself is expected to be
 an external input.
- Previously, some signs were inverted when writing the output, which is incovenient for
 debugging. Most importantly, this can lead to serious problems when coupling. All signs
 are now defined according to the surface of reference elipsoid z=0 and the center of the
 Earth z = -r_earth.
- The elastic solution is coupled to the viscous one in a clean way through the column
 anomalies. To decouple them, simply set `rho_litho = 0.0`.
- The compressibility correction factor for effective viscosity is now only computed once
 at the beginning.

## Minor changes

Formatting:
- indent = 4 spaces now much more consistent
- commas are usually followed by spaces for legibility
- max column length is 94 characters for better rendering in github... etc.
- including fftw was done repeatedly across many subroutines. It is now included at
 the beginning of each module using it.

Misc:
- Some names have been made more explicit. Example: `rho_a` is now `rho_uppermantle`.
- Any future task is tagged with "TODO" as comment. This combines well with the TODOS
 extension of VScode and allows spotting tasks with a simple search across documents,
 even without the VScode extension.
- Removed the option of static load, which is very specific to benchmark tests and
 makes the code longer / less intelegible. The performance advantage we were gaining
 from the static load option is minor now that the computation is accelerated.
- Most parameter files were previously in `output/`, which did not make sense.
 They are now gathered in `par/`.
- `README.md` slightly modified to include last functionalities.
- We now distinguish between diagnostics and prognostics to control the update frequency
 of the variables. Previously, this was controlled by `update_equil`, which has been
 renamed to `update_diagnostics` (almost always less frequent than prognostics).

Output:
- 2D fields are now output as 2D fields (and not 3D anymore).
- Reduced output to minimum for debugging.