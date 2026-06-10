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

  /* --- physics --- */
  config const gam   = UNSET_R;
  config const mu    = UNSET_R;
  config const grav1 = UNSET_R, grav2 = UNSET_R, grav3 = UNSET_R;
  config const rhoFloor = UNSET_R, prsFloor = UNSET_R;

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

  inline proc isUnset(v: int)    do return v == UNSET_I;
  inline proc isUnset(v: real)   do return isNan(v);
  inline proc isUnset(v: string) do return v == UNSET_S;
}
