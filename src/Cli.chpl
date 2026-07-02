/* Cli.chpl — raw command-line layer.
 *
 * Every runtime option is declared here as a `config const` whose
 * default is an "unset" sentinel.  Effective values are resolved in
 * Params.chpl with the precedence
 *
 *     command line  >  runtime_params.ini  >  built-in default
 *
 * so `--nx1=256` etc. behave exactly like ordinary Chapel config
 * consts, while anything not given on the command line may come from
 * the parameter file.
 */
module Cli {
  use Math;

  param UNSET_I = min(int);
  param UNSET_S = "\x01unset";
  const UNSET_R = nan;          // detected with isNan()

  /* parameter file location; the default file is optional, an
     explicitly given one must exist */
  config const paramsFile = UNSET_S;

  /* --- problem selection and grid --- */
  config const problem  = UNSET_S;
  config const geometry = UNSET_S;
  config const nx1 = UNSET_I, nx2 = UNSET_I, nx3 = UNSET_I;
  config const x1min = UNSET_R, x1max = UNSET_R;
  config const x2min = UNSET_R, x2max = UNSET_R;
  config const x3min = UNSET_R, x3max = UNSET_R;
  // grid law per direction: uniform | log | log-dec | stretch
  config const gridX1 = UNSET_S, gridX2 = UNSET_S, gridX3 = UNSET_S;
  config const stretchX1 = UNSET_R, stretchX2 = UNSET_R,
               stretchX3 = UNSET_R;          // geometric-progression ratio
  // cells of the uniform block anchoring a stretched direction
  config const stretchUniX1 = UNSET_I, stretchUniX2 = UNSET_I,
               stretchUniX3 = UNSET_I;

  /* --- physics --- */
  config const gam   = UNSET_R;
  config const eos   = UNSET_S;
  config const csIso = UNSET_R, csSlope = UNSET_R;
  config const mu    = UNSET_R;
  config const kappa = UNSET_R;
  config const grav1 = UNSET_R, grav2 = UNSET_R, grav3 = UNSET_R;
  config const gravCentral = UNSET_R, gravEps = UNSET_R;
  config const rhoFloor = UNSET_R, prsFloor = UNSET_R;
  // optically thin cooling: Lambda(T) = coolLambda0 * T^coolAlpha
  config const coolLambda0 = UNSET_R, coolAlpha = UNSET_R,
               coolTfloor = UNSET_R;
  config const scDiff = UNSET_R;        // passive-scalar diffusivity
  // spectral (Ornstein-Uhlenbeck) turbulence driving
  config const forceAmp = UNSET_R, forceTcorr = UNSET_R;
  config const forceKmin = UNSET_I, forceKmax = UNSET_I;
  config const forceSeed = UNSET_I;
  // Lagrangian tracer particles
  config const nParticles = UNSET_I, partSeed = UNSET_I;
  config const partRingR = UNSET_R;    // vortex ring seeding (example hook)
  // self-gravity: 4*pi*G normalisation (>0 enables the Poisson solve)
  config const sgFourPiG = UNSET_R, sgTol = UNSET_R;
  config const sgMaxIter = UNSET_I;
  // shearing box (Omega > 0 enables the rotating-frame source terms)
  config const omegaRot = UNSET_R, shearQ = UNSET_R;
  // FARGO orbital advection
  config const fargo = UNSET_S;

  /* --- numerics --- */
  config const cfl = UNSET_R, cflVisc = UNSET_R;
  config const tstop = UNSET_R, dtMax = UNSET_R;
  config const maxSteps = UNSET_I;
  config const recon = UNSET_S, limiter = UNSET_S;
  config const riemann = UNSET_S, integrator = UNSET_S;

  /* --- boundary conditions --- */
  config const bcX1min = UNSET_S, bcX1max = UNSET_S;
  config const bcX2min = UNSET_S, bcX2max = UNSET_S;
  config const bcX3min = UNSET_S, bcX3max = UNSET_S;

  /* --- output --- */
  config const outDt = UNSET_R;
  config const outFormats = UNSET_S, outDir = UNSET_S;
  config const logEvery = UNSET_I;
  /* resume from <outDir>/restart.chaa (CLI-only; see the docs'
     "Stopping & restarting" page) */
  config const restart = false;
  // periodic restart-dump cadence in simulation time (0 = off)
  config const restartDt = UNSET_R;

  /* --- problem-specific parameters --- */
  config const sodX0 = UNSET_R;
  config const sodRhoL = UNSET_R, sodVxL = UNSET_R, sodPrsL = UNSET_R;
  config const sodRhoR = UNSET_R, sodVxR = UNSET_R, sodPrsR = UNSET_R;
  config const cen1 = UNSET_R, cen2 = UNSET_R, cen3 = UNSET_R;
  config const sedovE0 = UNSET_R, sedovR0 = UNSET_R;
  config const sedovRhoAmb = UNSET_R, sedovPrsAmb = UNSET_R;
  config const blastPin = UNSET_R, blastPout = UNSET_R;
  config const blastRhoIn = UNSET_R, blastRhoOut = UNSET_R, blastR0 = UNSET_R;
  config const inRho = UNSET_R, inVx1 = UNSET_R, inVx2 = UNSET_R,
               inVx3 = UNSET_R, inPrs = UNSET_R;
  config const tcOmegaIn = UNSET_R, tcOmegaOut = UNSET_R;
  config const cylRad = UNSET_R;
  config const vortexBeta = UNSET_R;
  config const khRhoIn = UNSET_R, khRhoOut = UNSET_R, khV0 = UNSET_R,
               khPert = UNSET_R, khPrs = UNSET_R;
  config const rtRhoTop = UNSET_R, rtRhoBot = UNSET_R, rtPrs0 = UNSET_R,
               rtPert = UNSET_R;
  config const diskH0 = UNSET_R, diskJumpR = UNSET_R, diskJumpW = UNSET_R;
  config const twAmp = UNSET_R;
  config const cloudChi = UNSET_R, cloudRad = UNSET_R;
  config const waveAmp = UNSET_R;

  inline proc isUnset(v: int)    do return v == UNSET_I;
  inline proc isUnset(v: real)   do return isNan(v);
  inline proc isUnset(v: string) do return v == UNSET_S;
}
