# Tutorial

This walkthrough takes you from a 1D shock tube to a 3D blast wave in
spherical coordinates, ending in ParaView. Everything below assumes you
have [built Chaa](getting-started.md) (`build/bin/chaa`).

## 1. A shock tube, three ways

The Sod problem is the canonical first test. Run it with the default
HLLC solver and PLM reconstruction:

```sh
./build/bin/chaa --problem=sod --nx1=400 --tstop=0.2 --outFormats=txt
```

Now try the **exact (Godunov) Riemann solver** and **PPM**
reconstruction and compare:

```sh
./build/bin/chaa --problem=sod --nx1=400 --tstop=0.2 --outFormats=txt \
                 --riemann=exact --recon=ppm --integrator=rk3 \
                 --outDir=output-ppm
```

```python
import numpy as np, matplotlib.pyplot as plt
for d, lbl in (("output", "hllc + plm"), ("output-ppm", "exact + ppm")):
    x, rho = np.loadtxt(f"{d}/sod.0001.txt", usecols=(0, 1), unpack=True)
    plt.plot(x, rho, label=lbl)
plt.legend(); plt.show()
```

The same problem becomes a *radial* shock tube just by switching
geometry — note the reflecting axis boundary at r = 0:

```sh
./build/bin/chaa --problem=sod --geometry=spherical --nx1=400 \
                 --bcX1min=axis --tstop=0.2 --outFormats=txt
```

## 2. Sedov–Taylor in 2D (R,z)

A point explosion on an axisymmetric cylindrical mesh:

```sh
./build/bin/chaa --problem=sedov --geometry=cylindrical \
   --nx1=128 --nx2=256 --x1min=0 --x1max=1.2 --x2min=-1.2 --x2max=1.2 \
   --bcX1min=axis --sedovR0=0.04 --tstop=0.5 --outDt=0.1 \
   --outFormats=hdf5,vtk
```

This writes a time series `output/sedov.000N.h5` with `.xmf`
companions. The blast should be a perfect circle in the (R,z) plane —
the test suite checks it stays spherical to better than 0.1 %, with the
shock at the Sedov similarity radius
\( R(t) = (E t^2/\alpha\rho)^{1/5} \).

## 3. Into ParaView (or VisIt)

1. open `output/sedov.0005.xmf` (choose the *XDMF Reader* if asked),
2. press *Apply* — the mesh appears in physical coordinates (a half
   plane here; polar and spherical runs render as annuli, wedges and
   shells),
3. colour by `rho`, add a *Plot Over Line* along R for profiles.

The `.h5` files are also trivially scriptable:

```python
import h5py
with h5py.File("output/sedov.0005.h5") as f:
    rho = f["rho"][:]          # (nx2, nx1) — C order, x1 fastest
    R, z = f["cc_x1"][:], f["cc_x2"][:]
    t = f["time"][0]
```

## 4. A 3D spherical blast

```sh
./build/bin/chaa --problem=sedov --geometry=spherical \
   --nx1=64 --nx2=32 --nx3=32 \
   --x1min=0 --x1max=1.2 --x2min=0 --x2max=3.14159265 \
   --x3min=0 --x3max=6.28318530 \
   --bcX1min=reflect --bcX2min=axis --bcX2max=axis \
   --bcX3min=periodic --bcX3max=periodic \
   --sedovR0=0.12 --tstop=0.3 --outFormats=hdf5
```

In ParaView the 3DSMesh XDMF renders as a solid sphere of cells. The
shock radius is independent of θ and φ to machine precision — that is
the well-balanced geometry at work.

## 5. Physics switches in one line each

```sh
# isothermal EOS (sound speed 1), as in Idefix's sod-iso:
./build/bin/chaa --problem=sod --eos=isothermal --csIso=1 --nx1=500 --tstop=0.2

# viscous Taylor-Couette flow between rotating cylinders:
./build/bin/chaa --problem=taylorCouette --geometry=cylindrical --nx1=64 \
   --x1min=1 --x1max=2 --bcX1min=userdef --bcX1max=userdef \
   --mu=0.05 --inPrs=10 --tstop=30

# thermal conduction (decaying entropy mode):
./build/bin/chaa --problem=thermalWave --nx1=128 --x1min=-0.5 --x1max=0.5 \
   --bcX1min=periodic --bcX1max=periodic --kappa=0.02 --tstop=2

# Keplerian disk with a cavity (central gravity + locally isothermal):
./build/bin/chaa --problem=diskCavity --geometry=polar --nx1=96 --nx2=96 \
   --x1min=0.4 --x1max=2.5 --x2min=0 --x2max=6.2831853 \
   --bcX2min=periodic --bcX2max=periodic \
   --eos=isothermal --csIso=0.1 --csSlope=-0.5 --gravCentral=1 --tstop=10
```

## 6. Reproducible runs with a parameter file

Put your configuration in an INI file instead of a long command line:

```ini
# myrun.ini
[run]
problem  = blast
geometry = polar
tstop    = 0.15
[grid]
nx1 = 128
nx2 = 256
x1min = 0.5
x1max = 1.5
x2max = 6.28318530717959
[boundaries]
bcX2min = periodic
bcX2max = periodic
```

```sh
./build/bin/chaa --paramsFile=myrun.ini --cen1=-1 --outFormats=hdf5
```

Command-line flags always win over the file, which wins over built-in
defaults — see the [configuration guide](user-guide/configuration.md).

## Plot what you made

Every run in this tutorial can be inspected with the bundled tools:

```sh
python tools/plot_fields.py output/               # initial vs final fields
python tools/plot_compare.py sod output/          # overlay the exact solution
```

See [Plotting & analysis](user-guide/plotting.md) for all options.
