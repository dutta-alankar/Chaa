/* Chaa — a performance-portable finite-volume hydrodynamics solver
 * written in Chapel.
 *
 * Solves the compressible Euler / Navier-Stokes equations on uniform
 * structured grids in 1D/2D/3D Cartesian, cylindrical (axisymmetric),
 * polar and spherical coordinates.  All distribution and communication
 * is delegated to Chapel's StencilDist; the numerics below are written
 * as plain data-parallel foralls.
 *
 *   build :  make            (see Makefile; HDF5=1 enables HDF5 output)
 *   run   :  ./bin/chaa --problem=sod --nx1=400 --tstop=0.2 ...
 */
module Chaa {
  use Params, Grid, State, Eos, Hydro, Boundary, Problems, Evolve, Output;
  use Time, FileSystem, Math;

  proc printBanner() {
    writeln("=========================================================");
    writeln("  Chaa — Chapel hydrodynamics  (brewed fresh, served hot)");
    writeln("=========================================================");
    writeln("  problem    : ", problem);
    writeln("  geometry   : ", geometry, "  (", ndim, "D)");
    writeln("  grid       : ", nx1, " x ", nx2, " x ", nx3);
    writeln("  domain     : [", x1min, ",", x1max, "] x [",
            x2min, ",", x2max, "] x [", x3min, ",", x3max, "]");
    writeln("  solver     : ", riemann, " / ", recon, " (", limiter,
            ") / ", integrator);
    writeln("  gamma      : ", gam, "   cfl: ", cfl);
    if mu > 0.0 then writeln("  viscosity  : mu = ", mu);
    writeln("  locales    : ", numLocales, "   tasks/locale: ",
            here.maxTaskPar);
    writeln("=========================================================");
  }

  proc sanityChecks() {
    if !act1 then halt("nx1 must be > 1");
    if act3 && !act2 then halt("activate x2 before x3 (use nx2 > 1)");
    if geom == Geom.cylindrical && act3 then
      halt("cylindrical (R,z) is axisymmetric: nx3 must be 1; " +
           "use geometry=polar for (R,phi,z)");
    if mu > 0.0 && geom == Geom.spherical then
      halt("viscosity is not implemented in spherical coordinates");
    if geom == Geom.spherical && act2 &&
       (x2min < -1.0e-12 || x2max > pi + 1.0e-12) then
      halt("spherical: theta must lie in [0, pi]");
  }

  proc main() {
    printBanner();
    sanityChecks();
    if !(try! exists(outDir)) then try! mkdir(outDir, parents=true);

    problemInit();
    forall idx in DInt do U[idx] = prim2cons(V[idx]);
    problemInternalBC(0.0);
    applyBCs(0.0);

    var t = 0.0;
    var step = 0;
    var dumpN = 0;
    var nextOut = outDt;
    var sw: stopwatch;
    sw.start();

    writeOutputs(dumpN, t);
    dumpN += 1;

    var dt = computeDt();
    while t < tstop*(1.0 - 1.0e-12) && step < maxSteps {
      dt = min(dt, tstop - t);
      advance(dt, t);
      t += dt;
      step += 1;
      if step % logEvery == 0 then
        writeln("step ", step, "  t = ", t, "  dt = ", dt);
      dt = computeDt();
      if dt <= 0.0 || dt != dt then halt("invalid time step: ", dt);
      if outDt > 0.0 && t >= nextOut - 1.0e-12 &&
         t < tstop*(1.0 - 1.0e-12) {
        writeOutputs(dumpN, t);
        dumpN += 1;
        nextOut += outDt;
      }
    }

    writeOutputs(dumpN, t);
    sw.stop();

    const ncells = nx1*nx2*nx3;
    writeln("done: ", step, " steps to t = ", t, " in ",
            sw.elapsed(), " s  (",
            (step:real * ncells:real)/max(sw.elapsed(), 1e-30)/1.0e6,
            " Mcell-updates/s)");
  }
}
