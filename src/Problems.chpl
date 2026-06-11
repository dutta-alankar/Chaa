/* Problems.chpl — registry/dispatcher for the bundled test problems.
 *
 * Each problem lives in its own file under src/problems/ and provides
 *   setup()            initial condition (fills V over the interior)
 *   userBC(side, t)    optional, for --bcX*=userdef sides
 *   internalBC(t)      optional, immersed/internal boundaries
 *
 * To add a problem: create src/problems/MyProblem.chpl with a setup()
 * proc and register it in the three dispatchers below.
 */
module Problems {
  use Params;
  import Sod, TwoBlast, Sedov, Blast, Riemann2D, DoubleMach,
         KelvinHelmholtz, RayleighTaylor, IsentropicVortex,
         TaylorCouette, CylinderFlow, KhIdefix, DiskCavity, ThermalWave,
         Cloud, LinearWave, Turbulence, Epicycle;

  proc problemInit() {
    select problem {
      when "sod"           do Sod.setup();
      when "twoblast"      do TwoBlast.setup();
      when "sedov"         do Sedov.setup();
      when "blast"         do Blast.setup();
      when "riemann2d"     do Riemann2D.setup();
      when "dmr"           do DoubleMach.setup();
      when "kh"            do KelvinHelmholtz.setup();
      when "rt"            do RayleighTaylor.setup();
      when "vortex"        do IsentropicVortex.setup();
      when "taylorCouette" do TaylorCouette.setup();
      when "cylinderFlow"  do CylinderFlow.setup();
      when "khi"           do KhIdefix.setup();
      when "diskCavity"    do DiskCavity.setup();
      when "thermalWave"   do ThermalWave.setup();
      when "cloud"         do Cloud.setup();
      when "linearWave"    do LinearWave.setup();
      when "turbulence"    do Turbulence.setup();
      when "epicycle"      do Epicycle.setup();
      otherwise do halt("unknown problem: " + problem);
    }
  }

  proc problemUserBC(side: int, t: real) {
    select problem {
      when "dmr"           do DoubleMach.userBC(side, t);
      when "taylorCouette" do TaylorCouette.userBC(side, t);
      otherwise do halt("problem " + problem +
                        " does not define user boundary conditions");
    }
  }

  /* optional problem-defined body force (acceleration); register a
     problem here and return its acceleration vector — see the docs'
     custom-problem guide for a worked example */
  const problemHasBodyForce = false;     // || problem == "myproblem"

  proc problemBodyForce(i: int, j: int, k: int, t: real): 3*real {
    // select problem { when "myproblem" do return MyProblem.accel(i,j,k,t); }
    return (0.0, 0.0, 0.0);
  }

  proc problemInternalBC(t: real) {
    select problem {
      when "cylinderFlow" do CylinderFlow.internalBC(t);
      otherwise { }
    }
  }
}
