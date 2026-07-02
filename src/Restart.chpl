/* Restart.chpl — graceful stop and machine-identical restart.
 *
 * A restart file (<outDir>/restart.chaa, little-endian binary) captures
 * everything the time loop depends on:
 *   - the conservative state U over the interior,
 *   - t, step count, dump counter and next scheduled output time,
 *   - tracer-particle ids and positions,
 *   - the OU forcing amplitudes and the RNG draw count (so the random
 *     sequence resumes exactly where it left off),
 *   - the self-gravity potential (the warm start of the CG solver).
 *
 * With the same binary and launch configuration, a run resumed with
 * --restart=true therefore reproduces the uninterrupted run bit for
 * bit, including all subsequent dumps (validated by the `restart-sod`
 * CI case).  The file is written when a graceful stop is requested
 * (`touch <outDir>/stop`) and at every normal end of run.
 */
module Restart {
  use Params, State, Particles, Forcing, SelfGravity;
  use IO;

  param RST_MAGIC = 0x4348414152535431;   // "CHAARST1"

  proc restartPath(): string do return outDir + "/restart.chaa";

  proc writeRestart(t: real, step: int, dumpN: int,
                    nextOut: real) throws {
    {
      var f = open(restartPath(), ioMode.cw);
      var w = f.writer(locking=false);
      w.writeBinary(RST_MAGIC);
      w.writeBinary(NTOT); w.writeBinary(nx1);
      w.writeBinary(nx2); w.writeBinary(nx3);
      w.writeBinary(t); w.writeBinary(step);
      w.writeBinary(dumpN); w.writeBinary(nextOut);

      // conservative state, one (j,k) row at a time
      var buf: [0..#(nx1*NTOT)] real;
      for k in 1..nx3 do for j in 1..nx2 {
        forall i in 1..nx1 with (ref buf) do
          for param c in 0..NTOT-1 do buf[(i-1)*NTOT + c] = U[i,j,k](c);
        w.writeBinary(buf);
      }

      // tracer particles (id + position)
      w.writeBinary(nParticles);
      if nParticles > 0 {
        for l in 0..#numLocales {
          const arr = bag[l].toArray();
          for pt in arr {
            w.writeBinary(pt.id);
            for param c in 0..2 do w.writeBinary(pt.pos(c));
          }
        }
      }

      // OU forcing state
      w.writeBinary(if forceAmp > 0.0 then 1 else 0);
      if forceAmp > 0.0 {
        w.writeBinary(nModes);
        for m in 0..#nModes { w.writeBinary(am1[m]); w.writeBinary(am2[m]); }
        w.writeBinary(nDraws);
      }

      // self-gravity potential (warm start of the CG solver)
      w.writeBinary(if sgFourPiG > 0.0 then 1 else 0);
      if sgFourPiG > 0.0 {
        var gb: [0..#nx1] real;
        for k in 1..nx3 do for j in 1..nx2 {
          forall i in 1..nx1 with (ref gb) do gb[i-1] = PHI[i,j,k];
          w.writeBinary(gb);
        }
      }
      w.close();
    }
    writeln("restart state written to ", restartPath());
  }

  proc readRestart(ref t: real, ref step: int, ref dumpN: int,
                   ref nextOut: real) throws {
    {
      var f = open(restartPath(), ioMode.r);
      var r = f.reader(locking=false);
      var magic, ntot, n1, n2, n3: int;
      if !r.readBinary(magic) || magic != RST_MAGIC then
        halt("not a Chaa restart file: " + restartPath());
      r.readBinary(ntot); r.readBinary(n1); r.readBinary(n2);
      r.readBinary(n3);
      if ntot != NTOT || n1 != nx1 || n2 != nx2 || n3 != nx3 then
        halt("restart file was written for a ", n1, "x", n2, "x", n3,
             " grid with ", ntot, " fields; this run is ",
             nx1, "x", nx2, "x", nx3, " with ", NTOT);
      r.readBinary(t); r.readBinary(step);
      r.readBinary(dumpN); r.readBinary(nextOut);

      var buf: [0..#(nx1*NTOT)] real;
      for k in 1..nx3 do for j in 1..nx2 {
        r.readBinary(buf);
        forall i in 1..nx1 {
          var u: StateVec;
          for param c in 0..NTOT-1 do u(c) = buf[(i-1)*NTOT + c];
          U[i,j,k] = u;
        }
      }

      var np: int;
      r.readBinary(np);
      if np != nParticles then
        halt("restart file has ", np, " particles; this run asks for ",
             nParticles, " (--nParticles must match)");
      if np > 0 {
        var ids: [0..#np] int;
        var pps: [0..#np] 3*real;
        for p in 0..#np {
          r.readBinary(ids[p]);
          for param c in 0..2 do r.readBinary(pps[p](c));
        }
        restoreParticles(ids, pps);
      }

      var hasF: int;
      r.readBinary(hasF);
      if hasF != (if forceAmp > 0.0 then 1 else 0) then
        halt("restart file and run disagree on turbulence forcing");
      if hasF == 1 {
        var nm: int;
        r.readBinary(nm);
        if nm != nModes then
          halt("restart file has ", nm, " forcing modes, this run built ",
               nModes);
        var a1, a2: [0..#nm] real;
        for m in 0..#nm { r.readBinary(a1[m]); r.readBinary(a2[m]); }
        var draws: int;
        r.readBinary(draws);
        restoreForcing(a1, a2, draws);
      }

      var hasG: int;
      r.readBinary(hasG);
      if hasG != (if sgFourPiG > 0.0 then 1 else 0) then
        halt("restart file and run disagree on self-gravity");
      if hasG == 1 {
        var gb: [0..#nx1] real;
        for k in 1..nx3 do for j in 1..nx2 {
          r.readBinary(gb);
          forall i in 1..nx1 do PHI[i,j,k] = gb[i-1];
        }
      }
      r.close();
    }
    writeln("restarted from ", restartPath(), ": t = ", t, ", step ", step);
  }
}
