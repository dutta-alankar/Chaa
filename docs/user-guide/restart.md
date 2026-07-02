# Stopping & restarting

Any Chaa run can be stopped gracefully at any moment and resumed later
— and the resumed run reproduces the uninterrupted one **bit for bit**.

## Stopping a running simulation

Create an empty file named `stop` in the run's output directory:

```sh
touch <outDir>/stop        # e.g.  touch output/stop
```

At the end of the step it is currently computing, the simulation

1. writes a restart file, `<outDir>/restart.chaa`,
2. **removes the `stop` file**, and
3. exits with a log message telling you how to resume:

```
restart state written to output/restart.chaa
stop requested: state saved after step 1731 (t = 0.8422); resume with --restart=true
```

No signal handling or job-control tricks are involved — this works the
same on a laptop, inside a batch job, and in multi-locale runs.

A restart file is **also written at every normal end of run**, so a
finished simulation can be continued to a later `--tstop` without
planning ahead.

## Restarting

Re-run with the *same* flags plus `--restart=true`:

```sh
./build/bin/chaa --problem=turbulence ... --outDir=output              # stopped early
./build/bin/chaa --problem=turbulence ... --outDir=output --restart=true
```

The run resumes from the saved step: output files continue with the
same numbering and land at the same simulation times, so the directory
ends up exactly as if the run had never been interrupted.

## What the restart file contains

`restart.chaa` is a small (grid-sized) little-endian binary capturing
everything the time loop depends on:

| saved | why |
|---|---|
| conservative state `U` (interior) | the full field state |
| `t`, step count, dump counter, next output time | continuation + identical output cadence/numbering |
| tracer-particle ids and positions | particles resume their trajectories |
| OU forcing amplitudes + RNG draw count | the turbulence-driving random sequence continues *exactly where it left off* (the stream is fast-forwarded on restart) |
| self-gravity potential Φ | the CG warm start matches, so the iteration history is identical |

Grid size, field count, particle count and the enabled
forcing/self-gravity modules are checked on read; a mismatch aborts
with a clear message rather than silently mixing configurations.

## Machine-identical continuation

Provided Chaa **has not been recompiled** in between and the run is
resumed with the same launch configuration (same flags, locale count
and threads), the restarted run produces *identical* results — every
subsequent dump is byte-for-byte equal to the uninterrupted run's.
This is enforced in CI:

- `restart-sod` — a Sod run is stopped after its first step (a
  pre-planted `stop` file), resumed, and every dump is required to be
  **identical** to an uninterrupted reference run;
- `restart-turbulence` — the same protocol on the OU-driven turbulence
  box with tracer particles: the hardest case, since it also requires
  the random forcing sequence and the particle trajectories to
  continue exactly.

Caveats worth knowing:

- a recompile (different compiler version, flags, or optimisation
  level) may change floating-point instruction scheduling — the resumed
  run is then still *correct*, just not bit-identical;
- restarting is exact because the `stop` file interrupts *between*
  steps. Runs that ended at `--tstop` clamped their last `dt` to hit
  the stop time exactly, so continuing a *finished* run to a larger
  `--tstop` is a new trajectory after that point (correct, but not the
  same dt sequence as one long run would have taken);
- the restart file is overwritten in place on each save; copy it aside
  if you want to keep a specific checkpoint.
