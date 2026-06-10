/* Params.chpl — runtime configuration and global constants for Chaa.
 *
 * Every physical/numerical knob is a `config const`, settable on the
 * command line, e.g.  ./chaa --problem=sod --nx1=400 --tstop=0.2
 */
module Params {
  param NVAR = 5;

  // primitive variable slots
  param IRHO = 0, IVX1 = 1, IVX2 = 2, IVX3 = 3, IPRS = 4;
  // conservative variable slots
  param IMX1 = 1, IMX2 = 2, IMX3 = 3, IENG = 4;

  type StateVec = NVAR*real;

  enum Geom { cartesian, cylindrical, polar, spherical }

  /* --- problem selection and grid --- */
  config const problem  = "sod";
  config const geometry = "cartesian";   // cartesian | cylindrical | polar | spherical
  config const nx1 = 128, nx2 = 1, nx3 = 1;
  config const x1min = 0.0, x1max = 1.0;
  config const x2min = 0.0, x2max = 1.0;
  config const x3min = 0.0, x3max = 1.0;

  /* --- physics --- */
  config const gam      = 1.4;     // adiabatic index
  config const mu       = 0.0;     // dynamic viscosity coefficient
  config const grav1    = 0.0,     // constant gravity along coordinate axes
               grav2    = 0.0,
               grav3    = 0.0;
  config const rhoFloor = 1.0e-12, prsFloor = 1.0e-14;

  /* --- numerics --- */
  config const cfl        = 0.4;
  config const cflVisc    = 0.3;
  config const tstop      = 0.2;
  config const dtMax      = 1.0e30;
  config const maxSteps   = 1000000000;
  config const recon      = "linear";   // constant | linear
  config const limiter    = "vanleer";  // minmod | vanleer | mc
  config const riemann    = "hllc";     // llf | hll | hllc
  config const integrator = "rk2";      // euler | rk2 | rk3

  /* --- boundary conditions: periodic|outflow|reflect|axis|inflow|userdef --- */
  config const bcX1min = "outflow", bcX1max = "outflow";
  config const bcX2min = "outflow", bcX2max = "outflow";
  config const bcX3min = "outflow", bcX3max = "outflow";

  /* --- output --- */
  config const outDt      = 0.0;        // <=0 : dump only initial & final states
  config const outFormats = "txt";      // comma list among txt,vtk,hdf5
  config const outDir     = "output";
  config const logEvery   = 100;

  /* --- problem-specific parameters --- */
  // sod
  config const sodX0 = 0.5;
  config const sodRhoL = 1.0,   sodVxL = 0.0, sodPrsL = 1.0;
  config const sodRhoR = 0.125, sodVxR = 0.0, sodPrsR = 0.1;
  // sedov / blast common: explosion centre in *physical* (mapped) coordinates
  config const cen1 = 0.0, cen2 = 0.0, cen3 = 0.0;
  config const sedovE0 = 1.0, sedovR0 = 0.1;
  config const sedovRhoAmb = 1.0, sedovPrsAmb = 1.0e-5;
  // blast
  config const blastPin = 10.0, blastPout = 0.1;
  config const blastRhoIn = 1.0, blastRhoOut = 1.0, blastR0 = 0.1;
  // generic inflow state (inflow BC, flow past cylinder, dmr left state ...)
  config const inRho = 1.0, inVx1 = 0.0, inVx2 = 0.0, inVx3 = 0.0, inPrs = 1.0;
  // Taylor-Couette
  config const tcOmegaIn = 1.0, tcOmegaOut = 0.0;
  // flow past cylinder
  config const cylRad = 0.5;
  // isentropic vortex
  config const vortexBeta = 5.0;
  // Kelvin-Helmholtz
  config const khRhoIn = 2.0, khRhoOut = 1.0, khV0 = 0.5, khPert = 0.01, khPrs = 2.5;
  // Rayleigh-Taylor
  config const rtRhoTop = 2.0, rtRhoBot = 1.0, rtPrs0 = 2.5, rtPert = 0.01;

  /* --- derived --- */
  const geom = parseGeom(geometry);

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
  param RECON_CONST = 0, RECON_LINEAR = 1;
  param LIM_MINMOD = 0, LIM_VANLEER = 1, LIM_MC = 2;
  param RS_LLF = 0, RS_HLL = 1, RS_HLLC = 2;
  param TI_EULER = 0, TI_RK2 = 1, TI_RK3 = 2;
  param BC_PERIODIC = 0, BC_OUTFLOW = 1, BC_REFLECT = 2, BC_AXIS = 3,
        BC_INFLOW = 4, BC_USERDEF = 5;

  const reconCode = parseOpt(recon, ("constant", "linear"));
  const limCode   = parseOpt(limiter, ("minmod", "vanleer", "mc"));
  const rsCode    = parseOpt(riemann, ("llf", "hll", "hllc"));
  const tiCode    = parseOpt(integrator, ("euler", "rk2", "rk3"));

  const bcCode = (parseBC(bcX1min), parseBC(bcX1max),
                  parseBC(bcX2min), parseBC(bcX2max),
                  parseBC(bcX3min), parseBC(bcX3max));

  proc parseOpt(s: string, names): int {
    for param i in 0..<names.size do
      if s == names(i) then return i;
    halt("invalid option: " + s);
  }

  proc parseBC(s: string): int {
    select s {
      when "periodic" do return BC_PERIODIC;
      when "outflow"  do return BC_OUTFLOW;
      when "reflect"  do return BC_REFLECT;
      when "axis"     do return BC_AXIS;
      when "inflow"   do return BC_INFLOW;
      when "userdef"  do return BC_USERDEF;
      otherwise do halt("unknown boundary condition: " + s);
    }
  }
}
