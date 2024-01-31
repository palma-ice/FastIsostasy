# Refactoring for modularity, performance and sea-level treatment

This PR contains a substantial refactoring of the isostasy module. It does not comply
with the guideline of atomic PRs (one PR = one task) because of the time constraints.
However, all tests are passed by the refactor and a detailed description follows.

## Major changes

Modularity:
- `isostasy.f90` has been splitted into `convolutions.f90`, `finite_differences.f90`,
 `green_functions.f90`, `isos_utils.f90`, `kelvin_function.f90`, `sealevel.f90` and
 `solver_xlra.f90`. The latter includes the functions solving the viscous displacement
 of ELRA and LLRA (gathered under XLRA).
- Many functions were duplicated, which has been corrected
- The code that was used several times has now been written as functions
- `isostasy_defs.f90` defines a new derived type `isos_domain_class`, that stores all
 the information related to the computation domain (including the pseudo-differential
 operator in Fourier space as well as the FFT plans).


Performance and accuracy:
- FFT plans are precomputed
- Convolution kernels are precomputed
- Allocations inside the time loop are reduced
- The normalisation of the FFT is only applied on the backtransform, as it should be!

Sea-level:
- Column anomalies are now used as forcing

General physics:
- Allow rectangular domain (without the need of a square extension)
- Distortion factor included in computations


## Minor changes

Formatting:
- indent = 4 spaces now much more consistent
- commas are followed by spaces for legibility
- max column length is 94 characters for better rendering in github... etc.

Misc:
- Any future task is tagged with "TODO" as comment. This combines well with the TODOS
 extension of VScode and allows spotting tasks with a simple search across documents,
 even without the VScode extension.
- Removed the option of static load, which is very specific to benchmark tests and
 makes the code longer / less intelegible. The performance advantage we were gaining
 from the static load option is minor now that the computation is accelerated.