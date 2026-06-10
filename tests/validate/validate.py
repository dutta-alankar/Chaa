#!/usr/bin/env python3
"""Quantitative validation of Chaa test-problem output.

Usage:  validate.py <case> <outdir>

Each case checks the final dump of the corresponding run script in
tests/cases/ against exact solutions, similarity solutions, symmetry
requirements or documented reference values.
"""
import os
import sys
import glob
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import load_txt, load_h5, check, finish, sedov_radius
import exact_riemann


def last(outdir, ext):
    files = sorted(glob.glob(os.path.join(outdir, f"*.{ext}")))
    assert files, f"no .{ext} output in {outdir}"
    return files[-1]


# --------------------------------------------------------------------- #
def sod_1d_cart(outdir):
    d = load_txt(last(outdir, "txt"))
    x, rho, vx, prs = d[:, 0], d[:, 1], d[:, 2], d[:, 5]
    re, ue, pe = exact_riemann.solution(x, t=0.2)
    L1 = np.abs(rho - re).mean()
    print(f"  L1(rho) vs exact Riemann solution: {L1:.4e}")
    check(L1 < 0.012, f"L1(rho)={L1:.3e} < 0.012")
    check(np.abs(vx - ue).mean() < 0.025, "L1(vx) < 0.025")
    check(np.abs(prs - pe).mean() < 0.012, "L1(p) < 0.012")
    # hdf5 written for this case too: cross-check the two backends
    h = load_h5(last(outdir, "h5"))
    check(np.allclose(h["rho"], rho), "hdf5 rho matches txt rho")
    check(abs(h["time"][0] - 0.2) < 1e-10, "hdf5 time stamp == 0.2")


def sod_1d_exact_rs(outdir):
    """Godunov scheme with the built-in exact Riemann solver."""
    d = load_txt(last(outdir, "txt"))
    x, rho = d[:, 0], d[:, 1]
    re, ue, pe = exact_riemann.solution(x, t=0.2)
    L1 = np.abs(rho - re).mean()
    print(f"  L1(rho), exact-RS scheme vs exact solution: {L1:.4e}")
    check(L1 < 0.012, f"L1(rho)={L1:.3e} < 0.012")
    check(np.isfinite(d).all(), "all values finite")


def sod_1d_radial(outdir):
    d = load_txt(last(outdir, "txt"))
    rho = d[:, 1]
    check(np.isfinite(d).all(), "all values finite")
    check(rho.min() > 0.05 and rho.max() <= 1.0 + 1e-12,
          f"rho in (0.05, 1]: [{rho.min():.3f},{rho.max():.3f}]")
    # radial expansion weakens the jump relative to planar Sod
    jump = np.abs(np.diff(rho)).max()
    check(jump > 0.025, f"shock present (max jump {jump:.3f})")


def twoblast(outdir):
    d = load_txt(last(outdir, "txt"))
    x, rho = d[:, 0], d[:, 1]
    pk, xpk = rho.max(), x[rho.argmax()]
    print(f"  rho peak {pk:.3f} at x={xpk:.4f}")
    check(np.isfinite(d).all(), "all values finite")
    check(4.5 < pk < 7.5, f"peak density {pk:.2f} in (4.5, 7.5)")  # ref ~6
    check(0.70 < xpk < 0.82, f"peak position {xpk:.3f} in (0.70, 0.82)")


def sedov_1d_sph(outdir):
    d = load_txt(last(outdir, "txt"))
    r, rho = d[:, 0], d[:, 1]
    rad, ana = r[rho.argmax()], sedov_radius(0.5)
    err = abs(rad - ana) / ana
    print(f"  shock radius {rad:.4f}, analytic {ana:.4f}, err {err*100:.2f}%")
    check(err < 0.04, f"radius error {err*100:.2f}% < 4%")
    check(3.5 < rho.max() < 6.5, f"peak compression {rho.max():.2f}")


def sedov_2d_cyl(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, R, z = h["rho"], h["cc_x1"], h["cc_x2"]
    ana = sedov_radius(0.5)
    radR = R[rho[np.abs(z).argmin(), :].argmax()]
    radZ = abs(z[rho[:, 0].argmax()])
    print(f"  shock at R={radR:.4f} z={radZ:.4f}, analytic {ana:.4f}")
    check(abs(radR - ana) / ana < 0.05, "R-radius within 5%")
    check(abs(radR - radZ) / ana < 0.03, "spherical to 3% between R and z axes")
    check(np.isfinite(rho).all(), "all values finite")


def sedov_2d_sph(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, r = h["rho"], h["cc_x1"]
    radii = r[rho.argmax(axis=1)]
    ana = sedov_radius(float(h["time"][0]))
    err = abs(radii.mean() - ana) / ana
    print(f"  radius {radii.mean():.4f} +- {radii.std():.2e}, analytic {ana:.4f}")
    check(radii.std() / ana < 0.01, "shock radius independent of theta (<1%)")
    check(err < 0.05, f"radius error {err*100:.2f}% < 5%")


def sedov_3d_cart(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, x = h["rho"], h["cc_x1"]
    n = rho.shape[0] // 2
    rad = abs(x[rho[n, n, :].argmax()])
    ana = sedov_radius(0.5)
    err = abs(rad - ana) / ana
    print(f"  shock radius {rad:.3f}, analytic {ana:.4f}, err {err*100:.1f}%")
    check(err < 0.08, f"radius error {err*100:.1f}% < 8%")
    check(np.isfinite(rho).all(), "all values finite")
    # octant symmetry along the three axes
    lx = rho[n, n, :]
    ly = rho[n, :, n]
    lz = rho[:, n, n]
    check(np.allclose(lx, ly, rtol=1e-10) and np.allclose(lx, lz, rtol=1e-10),
          "axis profiles identical (x=y=z symmetry)")


def sedov_3d_sph(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, r = h["rho"], h["cc_x1"]
    radii = r[rho.argmax(axis=2)]
    ana = sedov_radius(float(h["time"][0]))
    print(f"  radius {radii.mean():.4f} spread {radii.std():.2e}, ana {ana:.4f}")
    check(np.isfinite(rho).all(), "all values finite")
    check(radii.std() / ana < 0.01, "radius independent of theta & phi")
    check(abs(radii.mean() - ana) / ana < 0.08, "radius within 8%")


def blast_2d_polar(outdir):
    h = load_h5(last(outdir, "h5"))
    rho = h["rho"]
    asym = np.abs(rho - rho[::-1, :]).max() / rho.max()
    print(f"  phi-mirror asymmetry {asym:.2e}")
    check(np.isfinite(rho).all(), "all values finite")
    check(asym < 1e-10, "mirror-symmetric about blast centre")
    check(rho.max() > 1.5 * rho.min(), "blast wave present")


def blast_3d_polar(outdir):
    h = load_h5(last(outdir, "h5"))
    rho = h["rho"]
    asym = np.abs(rho - rho[:, ::-1, :]).max() / rho.max()
    zsym = np.abs(rho - rho[::-1, :, :]).max() / rho.max()
    print(f"  phi asym {asym:.2e}, z asym {zsym:.2e}")
    check(np.isfinite(rho).all(), "all values finite")
    check(asym < 1e-10 and zsym < 1e-10, "phi and z mirror symmetry")
    # xmf companion must be valid XML
    import xml.etree.ElementTree as ET
    tree = ET.parse(last(outdir, "xmf"))
    check(tree.getroot().tag == "Xdmf", "xmf companion is valid XDMF XML")


def riemann2d(outdir):
    h = load_h5(last(outdir, "h5"))
    rho = h["rho"]
    dsym = np.abs(rho - rho.T).max() / rho.max()
    print(f"  diagonal asymmetry {dsym:.2e}, rho [{rho.min():.3f},{rho.max():.3f}]")
    check(np.isfinite(rho).all(), "all values finite")
    check(dsym < 1e-10, "symmetric across the x=y diagonal")
    check(0.05 < rho.min() and rho.max() < 2.2, "density within physical bounds")
    check(rho.max() > 1.6, "shock interaction developed (rho max > 1.6)")


def dmr(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, x = h["rho"], h["cc_x1"]
    row = rho[0, :]
    stem = x[np.where(row > 2.5)[0][-1]]
    print(f"  Mach stem foot at x={stem:.3f}, rho max {rho.max():.2f}")
    check(np.isfinite(rho).all(), "all values finite")
    check(2.4 < stem < 3.1, f"Mach stem foot {stem:.2f} in (2.4, 3.1)")
    check(rho.max() > 15.0, f"triple-point density {rho.max():.1f} > 15")


def kh(outdir):
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    g0 = np.abs(h0["vx2"]).max()
    g1 = np.abs(h["vx2"]).max()
    print(f"  max|vy|: {g0:.4f} -> {g1:.4f}")
    check(np.isfinite(h["rho"]).all(), "all values finite")
    check(g1 > 10 * g0, "shear instability grew by >10x")


def rt(outdir):
    h = load_h5(last(outdir, "h5"))
    rho, vy = h["rho"], h["vx2"]
    print(f"  max|vy| {np.abs(vy).max():.4f}")
    check(np.isfinite(rho).all(), "all values finite")
    check(np.abs(vy).max() > 0.02, "RT mode growing")
    check(rho.max() < 2.5 and rho.min() > 0.5, "densities bounded")


def vortex(outdir):
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    L1 = np.abs(h["rho"] - h0["rho"]).mean()
    drift = abs(h["rho"].mean() - h0["rho"].mean()) / h0["rho"].mean()
    print(f"  L1(rho) after one period: {L1:.3e}, mass drift {drift:.2e}")
    check(L1 < 8e-3, f"advection L1 error {L1:.2e} < 8e-3")
    check(drift < 1e-12, "mass conserved to machine precision")


def taylor_couette(outdir):
    d = load_txt(last(outdir, "txt"))
    R, vphi = d[:, 0], d[:, 4]
    R1, R2, Om1, Om2 = 1.0, 2.0, 1.0, 0.0
    a = (Om2 * R2**2 - Om1 * R1**2) / (R2**2 - R1**2)
    b = (Om1 - Om2) * R1**2 * R2**2 / (R2**2 - R1**2)
    ana = a * R + b / R
    err = np.abs(vphi - ana).mean() / np.abs(ana).mean()
    print(f"  L1 error vs analytic Couette profile: {err*100:.3f}%")
    check(err < 0.02, f"steady profile within 2% (got {err*100:.2f}%)")
    check(np.abs(d[:, 2]).max() < 0.1, "radial velocity small (steady state)")


def sod_1d_iso(outdir):
    """isothermal Sod tube vs the exact two-wave isothermal solution"""
    d = load_txt(last(outdir, "txt"))
    x, rho, vx, prs = d[:, 0], d[:, 1], d[:, 2], d[:, 5]
    re, ue = exact_riemann.iso_solution(x, t=0.2, cs=1.0)
    L1 = np.abs(rho - re).mean()
    print(f"  L1(rho) vs exact isothermal Riemann solution: {L1:.4e}")
    check(L1 < 0.012, f"L1(rho)={L1:.3e} < 0.012")
    check(np.abs(vx - ue).mean() < 0.03, "L1(vx) < 0.03")
    check(np.allclose(prs / rho, 1.0), "p = cs^2 rho enforced exactly")


def khi_iso(outdir):
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    g0 = np.abs(h0["vx2"]).max()
    g1 = np.abs(h["vx2"]).max()
    print(f"  max|vy|: {g0:.4f} -> {g1:.4f}")
    check(np.isfinite(h["rho"]).all(), "all values finite")
    check(g1 > 0.1, "shear instability grew from the interface seed")
    check(0.5 < h["rho"].min() and h["rho"].max() < 2.0,
          "subsonic isothermal: density stays near 1")


def thermal_diffusion(outdir):
    """entropy-mode decay rate vs Gamma = kappa (gam-1)/gam k^2"""
    d0 = load_txt(sorted(glob.glob(os.path.join(outdir, "*.txt")))[0])
    d1 = load_txt(last(outdir, "txt"))
    a0 = np.abs(d0[:, 1] - 1.0).max()
    a1 = np.abs(d1[:, 1] - 1.0).max()
    kap, gamma, kw, t = 0.02, 1.4, 2.0 * np.pi, 2.0
    ana = np.exp(-kap * (gamma - 1.0) / gamma * kw**2 * t)
    err = abs(a1 / a0 - ana) / ana
    print(f"  decay {a1/a0:.4f} vs analytic {ana:.4f} ({err*100:.2f}%)")
    check(err < 0.06, f"conduction decay rate within 6% (got {err*100:.2f}%)")


def disk_cavity(outdir):
    """locally isothermal disk + central gravity: equilibrium holds"""
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    R = h["cc_x1"]
    v0 = h0["vx2"].mean(axis=0)
    v1 = h["vx2"].mean(axis=0)
    sel = (R > 1.7) & (R < 2.2)          # outside the cavity jump
    drift = np.abs(v1[sel] - v0[sel]).max() / v0[sel].max()
    kep = np.abs(v1[sel] / R[sel] ** -0.5 - 1.0).max()
    print(f"  vphi drift {drift*100:.2f}%, deviation from Keplerian {kep*100:.2f}%")
    check(np.isfinite(h["rho"]).all(), "all values finite")
    check(drift < 0.02, "rotational equilibrium drift < 2%")
    check(kep < 0.03, "vphi within 3% of Keplerian outside the cavity")
    check(h["rho"].min() > 0.0, "density positive in the cavity")


def sedov_3d_idefix(outdir):
    """idefix SedovBlastWave: gamma=5/3 on [-0.5,0.5]^3, t=0.1"""
    h = load_h5(last(outdir, "h5"))
    rho, x = h["rho"], h["cc_x1"]
    n = rho.shape[0] // 2
    rad = abs(x[rho[n, n, :].argmax()])
    ana = sedov_radius(0.1, gamma=5.0 / 3.0)
    err = abs(rad - ana) / ana
    print(f"  shock radius {rad:.3f}, analytic {ana:.4f}, err {err*100:.1f}%")
    check(np.isfinite(rho).all(), "all values finite")
    check(err < 0.08, f"radius error {err*100:.1f}% < 8%")


def vortex_limo3(outdir):
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    L1 = np.abs(h["rho"] - h0["rho"]).mean()
    print(f"  L1(rho), limo3: {L1:.3e}  (linear reference ~2.1e-3)")
    check(L1 < 2.0e-3, f"limo3 beats linear reconstruction ({L1:.2e} < 2e-3)")


def vortex_ppm(outdir):
    h0 = load_h5(sorted(glob.glob(os.path.join(outdir, "*.h5")))[0])
    h = load_h5(last(outdir, "h5"))
    L1 = np.abs(h["rho"] - h0["rho"]).mean()
    print(f"  L1(rho), ppm+rk3: {L1:.3e}  (linear reference ~2.1e-3)")
    check(L1 < 6.0e-4, f"ppm is high-order accurate ({L1:.2e} < 6e-4)")


def cylinder_flow(outdir):
    h = load_h5(last(outdir, "h5"))
    vx, x, y = h["vx1"], h["cc_x1"], h["cc_x2"]
    X, Y = np.meshgrid(x, y)
    solid = X**2 + Y**2 < 0.5**2
    iy0 = np.abs(y).argmin()
    wake = vx[iy0, np.abs(x - 2.0).argmin()]
    print(f"  |v| inside solid: {np.abs(vx[solid]).max():.2e}, wake vx(2,0)={wake:.3f}")
    check(np.isfinite(vx).all(), "all values finite")
    check(np.abs(vx[solid]).max() < 1e-12, "no-slip solid enforced")
    check(wake < 0.15, f"wake deficit present (vx={wake:.3f} << 0.3)")
    # vtk written for this case: check legacy header
    with open(last(outdir, "vtk")) as f:
        head = [f.readline() for _ in range(5)]
    check(head[0].startswith("# vtk DataFile"), "legacy VTK header")
    check("RECTILINEAR_GRID" in "".join(head), "rectilinear VTK dataset")


CASES = {
    "sod-1d-cart": sod_1d_cart,
    "sod-1d-exact": sod_1d_exact_rs,
    "sod-1d-cyl": sod_1d_radial,
    "sod-1d-sph": sod_1d_radial,
    "twoblast-1d": twoblast,
    "sedov-1d-sph": sedov_1d_sph,
    "sedov-2d-cyl": sedov_2d_cyl,
    "sedov-2d-sph": sedov_2d_sph,
    "sedov-3d-cart": sedov_3d_cart,
    "sedov-3d-sph": sedov_3d_sph,
    "blast-2d-polar": blast_2d_polar,
    "blast-3d-polar": blast_3d_polar,
    "riemann2d": riemann2d,
    "dmr": dmr,
    "kh": kh,
    "rt": rt,
    "vortex": vortex,
    "taylor-couette": taylor_couette,
    "cylinder-flow": cylinder_flow,
    "sod-1d-iso": sod_1d_iso,
    "khi-2d-iso": khi_iso,
    "thermal-diffusion": thermal_diffusion,
    "disk-cavity": disk_cavity,
    "sedov-3d-idefix": sedov_3d_idefix,
    "vortex-limo3": vortex_limo3,
    "vortex-ppm": vortex_ppm,
}

if __name__ == "__main__":
    case, outdir = sys.argv[1], sys.argv[2]
    print(f"validating case '{case}' from {outdir}")
    CASES[case](outdir)
    finish()
