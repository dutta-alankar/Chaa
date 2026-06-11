/* Params.chpl — effective runtime parameters for Chaa.
 *
 * Each value is resolved with the precedence
 *     command line (--key=value, see Cli.chpl)
 *   > parameter file (runtime_params.ini, see IniReader.chpl)
 *   > built-in default (this file)
 * Compile-time parameters live in compile_params.chpl.
 */
module Params {
  import Cli;
  import IniReader;
  import CompileParams.NSCAL;

  param NVAR = 5;                 // hydro slots
  param NTOT = NVAR + NSCAL;      // hydro + passive tracer fields

  // primitive variable slots
  param IRHO = 0, IVX1 = 1, IVX2 = 2, IVX3 = 3, IPRS = 4;
  // conservative variable slots
  param IMX1 = 1, IMX2 = 2, IMX3 = 3, IENG = 4;
  // first passive-scalar slot (concentration in prims, rho*s in cons)
  param ISC = 5;

  type StateVec = NTOT*real;

  enum Geom { cartesian, cylindrical, polar, spherical }

  /* resolution helpers: command line > ini > default */
  proc valI(cli: int, key: string, def: int): int {
    if !Cli.isUnset(cli) then return cli;
    if IniReader.hasKey(key) then return IniReader.getI(key);
    return def;
  }
  proc valR(cli: real, key: string, def: real): real {
    if !Cli.isUnset(cli) then return cli;
    if IniReader.hasKey(key) then return IniReader.getR(key);
    return def;
  }
  proc valS(cli: string, key: string, def: string): string {
    if !Cli.isUnset(cli) then return cli;
    if IniReader.hasKey(key) then return IniReader.getS(key);
    return def;
  }

  /* --- problem selection and grid --- */
  const problem  = valS(Cli.problem,  "problem",  "sod");
  const geometry = valS(Cli.geometry, "geometry", "cartesian");
  const nx1 = valI(Cli.nx1, "nx1", 128),
        nx2 = valI(Cli.nx2, "nx2", 1),
        nx3 = valI(Cli.nx3, "nx3", 1);
  const x1min = valR(Cli.x1min, "x1min", 0.0),
        x1max = valR(Cli.x1max, "x1max", 1.0);
  const x2min = valR(Cli.x2min, "x2min", 0.0),
        x2max = valR(Cli.x2max, "x2max", 1.0);
  const x3min = valR(Cli.x3min, "x3min", 0.0),
        x3max = valR(Cli.x3max, "x3max", 1.0);
  // grid law per direction (uniform | log | log-dec | stretch)
  const gridX1 = valS(Cli.gridX1, "gridX1", "uniform"),
        gridX2 = valS(Cli.gridX2, "gridX2", "uniform"),
        gridX3 = valS(Cli.gridX3, "gridX3", "uniform");
  const stretchX1 = valR(Cli.stretchX1, "stretchX1", 1.05),
        stretchX2 = valR(Cli.stretchX2, "stretchX2", 1.05),
        stretchX3 = valR(Cli.stretchX3, "stretchX3", 1.05);

  /* --- physics --- */
  const gam      = valR(Cli.gam, "gam", 1.4);       // adiabatic index
  const eos      = valS(Cli.eos, "eos", "ideal");   // ideal | isothermal
  const csIso    = valR(Cli.csIso, "csIso", 1.0);   // isothermal sound speed
  const csSlope  = valR(Cli.csSlope, "csSlope", 0.0); // cs ~ csIso*R^csSlope
  const mu       = valR(Cli.mu, "mu", 0.0);         // dynamic viscosity
  const kappa    = valR(Cli.kappa, "kappa", 0.0);   // thermal conductivity
  const grav1    = valR(Cli.grav1, "grav1", 0.0),   // constant gravity
        grav2    = valR(Cli.grav2, "grav2", 0.0),
        grav3    = valR(Cli.grav3, "grav3", 0.0);
  const gravCentral = valR(Cli.gravCentral, "gravCentral", 0.0); // GM (point mass)
  const gravEps     = valR(Cli.gravEps, "gravEps", 0.0);         // softening
  const rhoFloor = valR(Cli.rhoFloor, "rhoFloor", 1.0e-12),
        prsFloor = valR(Cli.prsFloor, "prsFloor", 1.0e-14);
  // optically thin cooling: edot = -rho^2 Lambda0 T^alpha (Townsend exact)
  const coolLambda0 = valR(Cli.coolLambda0, "coolLambda0", 0.0),
        coolAlpha   = valR(Cli.coolAlpha, "coolAlpha", 0.5),
        coolTfloor  = valR(Cli.coolTfloor, "coolTfloor", 1.0e-6);
  const scDiff = valR(Cli.scDiff, "scDiff", 0.0);   // scalar diffusivity
  // turbulence driving (Ornstein-Uhlenbeck spectral forcing)
  const forceAmp   = valR(Cli.forceAmp, "forceAmp", 0.0),
        forceTcorr = valR(Cli.forceTcorr, "forceTcorr", 1.0);
  const forceKmin = valI(Cli.forceKmin, "forceKmin", 1),
        forceKmax = valI(Cli.forceKmax, "forceKmax", 2);
  const forceSeed = valI(Cli.forceSeed, "forceSeed", 1234);
  // Lagrangian tracer particles
  const nParticles = valI(Cli.nParticles, "nParticles", 0);
  const partSeed   = valI(Cli.partSeed, "partSeed", 4321);

  /* --- numerics --- */
  const cfl        = valR(Cli.cfl, "cfl", 0.4);
  const cflVisc    = valR(Cli.cflVisc, "cflVisc", 0.3);
  const tstop      = valR(Cli.tstop, "tstop", 0.2);
  const dtMax      = valR(Cli.dtMax, "dtMax", 1.0e30);
  const maxSteps   = valI(Cli.maxSteps, "maxSteps", 1000000000);
  const recon      = valS(Cli.recon, "recon", "linear");
  const limiter    = valS(Cli.limiter, "limiter", "vanleer");
  const riemann    = valS(Cli.riemann, "riemann", "hllc");
  const integrator = valS(Cli.integrator, "integrator", "rk2");

  /* --- boundary conditions: periodic | zero-gradient | reflect | axis |
         inflow | outflow-diode | inflow-diode | userdef --- */
  const bcX1min = valS(Cli.bcX1min, "bcX1min", "zero-gradient"),
        bcX1max = valS(Cli.bcX1max, "bcX1max", "zero-gradient");
  const bcX2min = valS(Cli.bcX2min, "bcX2min", "zero-gradient"),
        bcX2max = valS(Cli.bcX2max, "bcX2max", "zero-gradient");
  const bcX3min = valS(Cli.bcX3min, "bcX3min", "zero-gradient"),
        bcX3max = valS(Cli.bcX3max, "bcX3max", "zero-gradient");

  /* --- output --- */
  const outDt      = valR(Cli.outDt, "outDt", 0.0);
  const outFormats = valS(Cli.outFormats, "outFormats", "txt");
  const outDir     = valS(Cli.outDir, "outDir", "output");
  const logEvery   = valI(Cli.logEvery, "logEvery", 100);

  /* --- problem-specific parameters --- */
  const sodX0  = valR(Cli.sodX0, "sodX0", 0.5);
  const sodRhoL = valR(Cli.sodRhoL, "sodRhoL", 1.0),
        sodVxL  = valR(Cli.sodVxL, "sodVxL", 0.0),
        sodPrsL = valR(Cli.sodPrsL, "sodPrsL", 1.0);
  const sodRhoR = valR(Cli.sodRhoR, "sodRhoR", 0.125),
        sodVxR  = valR(Cli.sodVxR, "sodVxR", 0.0),
        sodPrsR = valR(Cli.sodPrsR, "sodPrsR", 0.1);
  const cen1 = valR(Cli.cen1, "cen1", 0.0),
        cen2 = valR(Cli.cen2, "cen2", 0.0),
        cen3 = valR(Cli.cen3, "cen3", 0.0);
  const sedovE0 = valR(Cli.sedovE0, "sedovE0", 1.0),
        sedovR0 = valR(Cli.sedovR0, "sedovR0", 0.1);
  const sedovRhoAmb = valR(Cli.sedovRhoAmb, "sedovRhoAmb", 1.0),
        sedovPrsAmb = valR(Cli.sedovPrsAmb, "sedovPrsAmb", 1.0e-5);
  const blastPin  = valR(Cli.blastPin, "blastPin", 10.0),
        blastPout = valR(Cli.blastPout, "blastPout", 0.1);
  const blastRhoIn  = valR(Cli.blastRhoIn, "blastRhoIn", 1.0),
        blastRhoOut = valR(Cli.blastRhoOut, "blastRhoOut", 1.0),
        blastR0     = valR(Cli.blastR0, "blastR0", 0.1);
  const inRho = valR(Cli.inRho, "inRho", 1.0),
        inVx1 = valR(Cli.inVx1, "inVx1", 0.0),
        inVx2 = valR(Cli.inVx2, "inVx2", 0.0),
        inVx3 = valR(Cli.inVx3, "inVx3", 0.0),
        inPrs = valR(Cli.inPrs, "inPrs", 1.0);
  const tcOmegaIn  = valR(Cli.tcOmegaIn, "tcOmegaIn", 1.0),
        tcOmegaOut = valR(Cli.tcOmegaOut, "tcOmegaOut", 0.0);
  const cylRad = valR(Cli.cylRad, "cylRad", 0.5);
  const vortexBeta = valR(Cli.vortexBeta, "vortexBeta", 5.0);
  const khRhoIn  = valR(Cli.khRhoIn, "khRhoIn", 2.0),
        khRhoOut = valR(Cli.khRhoOut, "khRhoOut", 1.0),
        khV0     = valR(Cli.khV0, "khV0", 0.5),
        khPert   = valR(Cli.khPert, "khPert", 0.01),
        khPrs    = valR(Cli.khPrs, "khPrs", 2.5);
  const rtRhoTop = valR(Cli.rtRhoTop, "rtRhoTop", 2.0),
        rtRhoBot = valR(Cli.rtRhoBot, "rtRhoBot", 1.0),
        rtPrs0   = valR(Cli.rtPrs0, "rtPrs0", 2.5),
        rtPert   = valR(Cli.rtPert, "rtPert", 0.01);
  const diskH0    = valR(Cli.diskH0, "diskH0", 0.1),     // disk aspect ratio
        diskJumpR = valR(Cli.diskJumpR, "diskJumpR", 1.5),
        diskJumpW = valR(Cli.diskJumpW, "diskJumpW", 0.15);
  const twAmp = valR(Cli.twAmp, "twAmp", 1.0e-3);   // thermalWave amplitude
  const cloudChi = valR(Cli.cloudChi, "cloudChi", 10.0),  // density contrast
        cloudRad = valR(Cli.cloudRad, "cloudRad", 0.25);
  const waveAmp = valR(Cli.waveAmp, "waveAmp", 1.0e-4);   // linearWave

  /* --- derived --- */
  const geom = parseGeom(geometry);

  /* output/dataset name of state-vector slot c */
  proc fieldName(c: int): string {
    select c {
      when 0 do return "rho";
      when 1 do return "vx1";
      when 2 do return "vx2";
      when 3 do return "vx3";
      when 4 do return "prs";
      otherwise do return "sc" + (c - ISC):string;
    }
  }

  proc parseGeom(s: string): Geom {
    select s {
      when "cartesian"   do return Geom.cartesian;
      when "cylindrical" do return Geom.cylindrical;
      when "polar"       do return Geom.polar;
      when "spherical"   do return Geom.spherical;
      otherwise do halt("unknown geometry: " + s);
    }
  }

  // which dimensions actually carry cells
  const act1 = nx1 > 1,
        act2 = nx2 > 1,
        act3 = nx3 > 1;
  const ndim = (if act1 then 1 else 0) + (if act2 then 1 else 0)
             + (if act3 then 1 else 0);

  // integer codes (avoid string compares in hot loops)
  param RECON_CONST = 0, RECON_LINEAR = 1, RECON_LIMO3 = 2, RECON_PPM = 3,
        RECON_WENOZ = 4;
  param LIM_MINMOD = 0, LIM_VANLEER = 1, LIM_MC = 2;
  param RS_LLF = 0, RS_HLL = 1, RS_HLLC = 2, RS_EXACT = 3;
  param TI_EULER = 0, TI_RK2 = 1, TI_RK3 = 2, TI_VL2 = 3;
  param EOS_IDEAL = 0, EOS_ISO = 1;
  param GRID_UNIFORM = 0, GRID_LOG = 1, GRID_LOGDEC = 2, GRID_STRETCH = 3;
  param BC_PERIODIC = 0, BC_ZEROGRAD = 1, BC_REFLECT = 2, BC_AXIS = 3,
        BC_INFLOW = 4, BC_USERDEF = 5, BC_OUT_DIODE = 6, BC_IN_DIODE = 7;

  const eosCode = parseOpt(eos, ("ideal", "isothermal"));
  const reconCode = parseOpt(recon,
                             ("constant", "linear", "limo3", "ppm", "wenoz"));
  const limCode   = parseOpt(limiter, ("minmod", "vanleer", "mc"));
  const rsCode    = parseOpt(riemann, ("llf", "hll", "hllc", "exact"));
  const tiCode    = parseOpt(integrator, ("euler", "rk2", "rk3", "vl2"));

  const gridCode = (parseGrid(gridX1), parseGrid(gridX2), parseGrid(gridX3));

  const bcCode = (parseBC(bcX1min), parseBC(bcX1max),
                  parseBC(bcX2min), parseBC(bcX2max),
                  parseBC(bcX3min), parseBC(bcX3max));

  proc parseOpt(s: string, names): int {
    for param i in 0..<names.size do
      if s == names(i) then return i;
    halt("invalid option: " + s);
  }

  proc parseGrid(s: string): int {
    select s {
      when "uniform"  do return GRID_UNIFORM;
      when "log"      do return GRID_LOG;      // spacing grows with x
      when "log-dec"  do return GRID_LOGDEC;   // spacing shrinks with x
      when "stretch"  do return GRID_STRETCH;  // geometric progression
      otherwise do halt("unknown grid law: " + s);
    }
  }

  proc parseBC(s: string): int {
    select s {
      when "periodic"      do return BC_PERIODIC;
      when "zero-gradient" do return BC_ZEROGRAD;
      when "outflow"       do return BC_ZEROGRAD;  // legacy alias
      when "reflect"       do return BC_REFLECT;
      when "axis"          do return BC_AXIS;
      when "inflow"        do return BC_INFLOW;
      when "outflow-diode" do return BC_OUT_DIODE; // zero-grad, no inflow
      when "inflow-diode"  do return BC_IN_DIODE;  // zero-grad, no outflow
      when "userdef"       do return BC_USERDEF;
      otherwise do halt("unknown boundary condition: " + s);
    }
  }
}
