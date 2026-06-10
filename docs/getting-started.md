# Getting started

## Prerequisites

| dependency | purpose | install |
|---|---|---|
| [Chapel](https://chapel-lang.org) ≥ 2.8 | the compiler | `brew install chapel` / [docs](https://chapel-lang.org/download.html) / `docker pull chapel/chapel:2.8.0` |
| CMake ≥ 3.16 | build orchestration | `brew install cmake` / `apt install cmake` |
| libhdf5 *(optional)* | HDF5 output | `brew install hdf5` / `apt install libhdf5-dev` |
| python3 + numpy + h5py *(optional)* | test-suite validation | `pip install numpy h5py` |

## Build

```sh
git clone https://github.com/dutta-alankar/Chaa.git
cd Chaa
cmake -B build
cmake --build build        # -> build/bin/chaa
```

Useful configure options:

```sh
cmake -B build -DCHAA_HDF5=OFF          # drop the HDF5 dependency
cmake -B build -DHDF5_ROOT=/opt/hdf5    # point at a specific HDF5
cmake -B build -DCHAA_NG=3              # ghost layers (3 = default, needed for ppm)
cmake -B build -DCHAA_CHPL_FLAGS="--fast --detailed-errors"
```

## First run

```sh
./build/bin/chaa --problem=sod --nx1=400 --tstop=0.2 --outFormats=txt,vtk
```

You will see the banner, periodic log lines, and dumps in `output/`:

```
output/sod.0000.txt   # initial state
output/sod.0001.txt   # final state (t = 0.2)
output/sod.0000.vtk   # same, as legacy VTK
...
```

The txt format is plain columns (`x1 rho vx1 vx2 vx3 prs`) — plot it with
anything:

```python
import numpy as np, matplotlib.pyplot as plt
x, rho = np.loadtxt("output/sod.0001.txt", usecols=(0, 1), unpack=True)
plt.plot(x, rho); plt.xlabel("x"); plt.ylabel(r"$\rho$"); plt.show()
```

## Run the test suite

```sh
ctest --test-dir build -j 4        # all 26 validated cases
tests/run_case.sh sedov-2d-cyl     # ... or any single case
```

Each case runs the solver at CI resolution and validates the result
quantitatively (exact Riemann solutions, Sedov similarity scaling,
analytic Couette profiles, conduction decay rates, …).

## Where next

- the [tutorial](tutorial.md) walks through 1D/2D/3D runs and ParaView,
- the [configuration guide](user-guide/configuration.md) lists every
  parameter,
- [set up your own problem](custom-problem.md) shows how to add physics
  setups of your own.
