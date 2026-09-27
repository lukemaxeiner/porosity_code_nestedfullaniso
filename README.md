# FiberForm Porosity And Covariance Tool

This folder contains a C program that reads binarized TIFF stacks with
`libtiff` and supports a staged workflow for porosity analysis, empirical
covariance construction, covariance-model fitting, and porosity-realization
generation.

## Current Capabilities

The code can perform the following tasks:

1. Compute global porosity from a binarized TIFF stack.
2. Compute a blockwise local porosity field using a shared
   `FILTER_RESOLUTION`.
3. Compute empirical porosity covariance from that porosity field.
4. Fit an isotropic exponential covariance model.
5. Generate continuous Gaussian porosity realizations on a user-specified
   simulation grid.

The code can run those stages individually or as a full pipeline.

## Input Convention

Each input stack is given as:

```txt
INPUT_DATA = "/path/to/stack.tif 2.0" ;
FILTER_RESOLUTION = 50.0
```

where:

- `INPUT_DATA` contains the TIFF path and that stack's `VOXEL_RESOLUTION`
- `FILTER_RESOLUTION` is shared across the run

Multiple stacks are supported:

```txt
INPUT_DATA = "/path/to/stack1.tif 2.0" ; "/path/to/stack2.tif 1.5" ;
FILTER_RESOLUTION = 50.0
```

## Run Modes

The current stage switches are:

```txt
RUN_FULL_PIPELINE = YES|NO
EMPIRICAL_COVARIANCE = YES|NO
FIT_COVARIANCE_MODEL = YES|NO
GENERATE_POROSITY_REALIZATIONS = YES|NO
```

Additional optional utilities still exist:

```txt
GLOBAL_POROSITY = YES|NO
POROSITY_FIELD = YES|NO
POROSITY_COVARIANCE_LAG = NO|3DLAG|CURVE
POROSITY_AXIS = X|Y|Z
POROSITY_COVARIANCE_REF_WINDOW = YES|NO
REF_WINDOW_VALUES = <x> <y> <z>
POROSITY_REALIZATION_COUNT = <N>
POROSITY_REALIZATION_GRID_RES = <Nx> <Ny> <Nz>
POROSITY_REALIZATION_GRID_DIM = <x> <y> <z>
```

## Output Locations

Empirical covariance and fit artifacts are written next to the source TIFF in:

```txt
BUILD_COVARIANCE
```

This includes items such as:

- `*_porosity_covariance_lag.txt`
- `*_porosity_covariance_3dlag_filter_<...>.vtk`
- `*_tail_noise_analysis.txt`
- `*_exponential_fit_quality.txt`
- `*_covariance_build_summary.txt`
- `*_fitted_covariance_model.txt`

Generated porosity realizations are written next to the source TIFF in:

```txt
GENERATED_POROSITY_FIELDS
```

Each realization is written as:

- floating-point TIFF stack
- legacy VTK volume for ParaView

## Porosity Convention

The active code path currently treats binary value `0` as void in the porosity
counting and porosity-field construction logic.

This is important:

- the live code is currently using `0 = void`
- if your binarized TIFF stacks use the opposite convention, the porosity logic
  must be adjusted

## Build

```bash
make
```

## Notes

- The realization generator currently uses a full covariance matrix and
  Cholesky factorization, so it is intended for modest simulation grids.
- The covariance fit is currently isotropic exponential.
- The local porosity field is currently computed on non-overlapping blocks, not
  a sliding window.
