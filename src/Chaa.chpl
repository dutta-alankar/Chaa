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
  use Forcing, Particles, SelfGravity, Fargo, Restart, Gpu;
  use Logo;
  use Time, FileSystem, Math;
  use CommDiagnostics, GpuDiagnostics;
  import Cli;
  import CompileParams.gpuEnabled;

  /* --commDiag=true: count remote gets/puts/on-stmts over the time
     loop and print them per locale at the end (multi-locale tuning) */
  config const commDiag = false;

  /* --gpuDiag=true (GPU build): count kernel launches and host<->device
     transfers over the time loop and print them at the end */
  config const gpuDiag = false;

  proc printBanner() {
    printLogo();
    writeln("=============================================================");
    writeln("  Chaa — Chapel-based Hydrodynamics for Astrophysical");
    writeln("         Applications            (চা: brewed fresh, served hot)");
    writeln("=============================================================");
    writeln("  problem    : ", problem);
    writeln("  geometry   : ", geometry, "  (", ndim, "D)");
    writeln("  grid       : ", nx1, " x ", nx2, " x ", nx3,
            "   (", gridX1, ", ", gridX2, ", ", gridX3, ")");
    writeln("  domain     : [", x1min, ",", x1max, "] x [",
            x2min, ",", x2max, "] x [", x3min, ",", x3max, "]");
    if NTOT > NVAR then writeln("  tracers    : ", NTOT - NVAR,
                                " passive scalar field(s)");
    if nParticles > 0 then writeln("  particles  : ", nParticles,
                                   " Lagrangian tracers");
    writeln("  solver     : ", riemann, " / ", recon, " (", limiter,
            ") / ", integrator);
    if eosCode == EOS_ISO then
      writeln("  eos        : isothermal (csIso=", csIso,
              ", csSlope=", csSlope, ")");
    else
      writeln("  eos        : ideal, gamma = ", gam);
    writeln("  cfl        : ", cfl);
    if mu > 0.0 then writeln("  viscosity  : mu = ", mu);
    if kappa > 0.0 then writeln("  conduction : kappa = ", kappa);
    if scDiff > 0.0 then writeln("  scalar diff: D = ", scDiff);
    if gravCentral > 0.0 then writeln("  gravity    : GM = ", gravCentral);
    if coolLambda0 > 0.0 then
      writeln("  cooling    : Lambda0 = ", coolLambda0,
              ", alpha = ", coolAlpha);
    if sgFourPiG > 0.0 then
      writeln("  self-grav  : 4 pi G = ", sgFourPiG);
    if omegaRot > 0.0 then
      writeln("  shear box  : Omega = ", omegaRot, ", q = ", shearQ);
    if useFargo then writeln("  fargo      : orbital advection on");
    writeln("  locales    : ", numLocales, "   tasks/locale: ",
            here.maxTaskPar);
    if gpuEnabled then
      writeln("  gpus       : ",
              + reduce ([loc in Locales] loc.gpus.size),
              " visible across ", numLocales, " locale(s)");
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
    if (reconCode == RECON_PPM || reconCode == RECON_WENOZ) && ng1 < 3 then
      halt("ppm/wenoz reconstruction needs 3 ghost layers: " +
           "rebuild with -DCHAA_NG=3 (the default)");
    if kappa > 0.0 && eosCode == EOS_ISO then
      halt("thermal conduction requires the ideal-gas EOS");
    if coolLambda0 > 0.0 && eosCode == EOS_ISO then
      halt("cooling requires the ideal-gas EOS");
    if (gridCode(0) == GRID_LOG || gridCode(0) == GRID_LOGDEC) &&
       x1min <= 0.0 then
      halt("log grids in x1 need x1min > 0");
    if (gridCode(1) == GRID_LOG || gridCode(1) == GRID_LOGDEC) &&
       x2min <= 0.0 then
      halt("log grids in x2 need x2min > 0");
    if (gridCode(2) == GRID_LOG || gridCode(2) == GRID_LOGDEC) &&
       x3min <= 0.0 then
      halt("log grids in x3 need x3min > 0");
    if forceAmp > 0.0 && geom != Geom.cartesian then
      halt("turbulence forcing is implemented for Cartesian boxes");
    if (bcCode(0) == BC_SHEAR || bcCode(1) == BC_SHEAR) &&
       (geom != Geom.cartesian || omegaRot <= 0.0) then
      halt("shear-periodic boundaries need a Cartesian shearing box " +
           "(--omegaRot > 0)");
    if useFargo {
      if !act2 then halt("fargo needs an active x2 direction");
      if gridCode(1) != GRID_UNIFORM then
        halt("fargo needs a uniform x2 grid");
      if geom == Geom.cartesian && omegaRot <= 0.0 then
        halt("fargo in Cartesian needs the shearing box (--omegaRot>0)");
      if geom == Geom.polar && gravCentral <= 0.0 then
        halt("fargo in polar needs --gravCentral > 0 (Keplerian w)");
      if geom == Geom.cylindrical || geom == Geom.spherical then
        halt("fargo is implemented for polar and shearing-box runs");
    }
    if geom == Geom.spherical && act2 &&
       (x2min < -1.0e-12 || x2max > pi + 1.0e-12) then
      halt("spherical: theta must lie in [0, pi]");
    if gpuEnabled {
      if problem == "cylinderFlow" then
        halt("cylinderFlow re-imposes internal (immersed) boundaries " +
             "every stage on the host — not supported in the GPU " +
             "build yet; use the CPU build");
      if problemHasBodyForce then
        halt("the problemBodyForce hook runs on the host — not " +
             "supported in the GPU build yet; use the CPU build");
    }
  }

  proc main(args: [] string) {
    Restart.exePath = args[0];
    printBanner();
    sanityChecks();
    if !(try! exists(outDir)) then try! mkdir(outDir, parents=true);

    problemInit();
    if eosCode == EOS_ISO then
      forall (i, j, k) in DInt do
        V[i,j,k](IPRS) = V[i,j,k](IRHO)*cs2At(i, j, k);
    forall idx in DInt do U[idx] = prim2cons(V[idx]);
    problemInternalBC(0.0);
    applyBCs(0.0);
    initForcing();
    initParticles();

    var t = 0.0;
    var step = 0;
    var dumpN = 0;
    var nextOut = outDt;

    if Cli.restart {
      try! readRestart(t, step, dumpN, nextOut);
      applyFloorsAndPrims();
      problemInternalBC(t);
      applyBCs(t);
    } else {
      writeOutputs(dumpN, t);
      dumpN += 1;
    }
    // periodic restart dumps: next multiple of restartDt after t
    var nextRst = if restartDt > 0.0
                  then restartDt*(floor(t/restartDt) + 1.0)
                  else 0.0;

    /* GPU build: carve the grid into per-device blocks and upload the
       (possibly restart-restored) state; from here on the devices own
       the interior and the host arrays are refreshed on demand */
    gpuInit();
    if gpuEnabled {
      /* every hot loop must actually run as a device kernel — an
         ineligible loop silently falls back to the host and reads
         device memory element-wise (~1000x slower, not wrong).  Count
         the launches of one full RHS + dt evaluation: per block
         that is 1 (zero) + 3 per active dim (recon+solve+div)
         + 1 (sources) + 2 (dt kernel + reduction). */
      startGpuDiagnostics();
      computeRHS(0.0);
      computeDt();
      stopGpuDiagnostics();
      const gd = getGpuDiagnostics();
      var launches = 0;
      for x in gd do launches += x.kernel_launch: int;
      const expected = (4 + 3*ndim)*nGpuTotal;
      if launches < expected then
        halt("GPU eligibility check failed: only ", launches,
             " kernel launches for one RHS+dt evaluation (expected ",
             expected, ") — a hot loop fell back to host execution; ",
             "recompile with --report-gpu to identify it");
      writeln("  gpu        : eligibility check ok (", launches,
              " kernel launches per RHS+dt across ", nGpuTotal,
              " device(s))");
    }

    var stopped = false;
    if commDiag then startCommDiagnostics();
    if gpuEnabled && gpuDiag then startGpuDiagnostics();
    /* the reported Mcell-updates/s measures the evolution loop only
       (initial/final output and restart I/O excluded) */
    var sw: stopwatch;
    sw.start();
    var dt = computeDt();
    while t < tstop*(1.0 - 1.0e-12) && step < maxSteps {
      dt = min(dt, tstop - t);
      updateForcing(dt);
      if gpuEnabled && forceAmp > 0.0 then gpuUpForcing();
      if gpuEnabled && sgFourPiG > 0.0 {
        // the CG solve stays on the host: pull the density, push the
        // refreshed potential to the devices
        gpuDownV();
        solveGravity();
        gpuUpPhi();
      } else {
        solveGravity();
      }
      advance(dt, t);
      if useFargo {
        // FARGO's row remap stages through the host in the GPU build
        if gpuEnabled then gpuDownU();
        fargoShift(dt);
        applyFloorsAndPrims();
        applyBCs(t + dt);
        if gpuEnabled then gpuUpAll();
      }
      if gpuEnabled && nParticles > 0 then gpuDownV();
      advanceParticles(dt, t);
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
      if restartDt > 0.0 && t >= nextRst - 1.0e-12 &&
         t < tstop*(1.0 - 1.0e-12) {
        try! writeRestart(t, step, dumpN, nextOut);
        nextRst += restartDt;
      }
      /* graceful stop: `touch <outDir>/stop` finishes this step, saves
         a restart file, removes the stop file and exits */
      if (try! exists(outDir + "/stop")) {
        try! writeRestart(t, step, dumpN, nextOut);
        try! remove(outDir + "/stop");
        writeln("stop requested: state saved after step ", step,
                " (t = ", t, "); resume with --restart=true");
        stopped = true;
        break;
      }
    }

    sw.stop();
    if commDiag {
      stopCommDiagnostics();
      const d = getCommDiagnostics();
      for l in 0..#numLocales do
        writeln("commDiag locale ", l, ": ", d[l]);
    }
    if gpuEnabled && gpuDiag {
      stopGpuDiagnostics();
      const d = getGpuDiagnostics();
      for l in 0..#numLocales do
        writeln("gpuDiag locale ", l, ": ", d[l]);
    }
    if gpuEnabled then gpuPrintTimers();
    if !stopped {
      writeOutputs(dumpN, t);
      // a failed restart write (quota, stale filesystem metadata)
      // must not kill an otherwise complete run
      try {
        writeRestart(t, step, dumpN + 1, nextOut);
      } catch e {
        writeln("WARNING: final restart write failed: ", e.message());
      }
    }

    const ncells = nx1*nx2*nx3;
    writeln("done: ", step, " steps to t = ", t, " in ",
            sw.elapsed(), " s  (",
            (step:real * ncells:real)/max(sw.elapsed(), 1e-30)/1.0e6,
            " Mcell-updates/s)");
  }
}
