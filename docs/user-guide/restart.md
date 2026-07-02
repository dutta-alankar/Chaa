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

## When restart files are written

A restart file is written

- on a **graceful stop** (the `stop` file, above),
- at every **normal end of run** — so a finished simulation can be
  continued to a later `--tstop` without planning ahead, and
- **periodically**, every `restartDt` of simulation time, if given in
  the parameter file or on the command line
  (`restartDt = 0.5` / `--restartDt=0.5`; default 0 = off) — insurance
  against crashes and killed batch jobs.

Writing a restart file never perturbs the run: a run that wrote
periodic restarts produces exactly the same dumps as one that didn't.

**Two generations are kept.** Before a new restart file is written, an
existing `restart.chaa` is renamed to `restart.bak.chaa` (replacing any
older backup). If the newest restart is unusable — a job killed
mid-write, say — rename the `.bak` file back and resume from one
checkpoint earlier.

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
| fingerprints of the **Chaa binary** and the **parameter file** | change detection at restart (below) |

Grid size, field count, particle count and the enabled
forcing/self-gravity modules are checked on read; a mismatch aborts
with a clear message rather than silently mixing configurations.

## Change detection

Every restart file records FNV-1a fingerprints of the Chaa executable
and of the parameter file in use when it was written. On
`--restart=true` the current fingerprints are compared and the result
is logged:

```
restart fingerprints match: same Chaa binary and parameter file — continuation is machine-identical
```

or, if something changed in between:

```
WARNING: the Chaa binary differs from the one that wrote this restart file — ...
WARNING: the parameter file differs from the one used when this restart file was written — ...
```

A mismatch is a warning, not an error — resuming with a retuned
parameter file or a rebuilt binary is a legitimate thing to do — but
the bit-identity guarantee below only holds when both fingerprints
match.

## Machine-identical continuation

Provided the **same Chaa binary and parameter file** are used (the
fingerprints confirm this in the log) and the run is resumed with the
same launch configuration (same flags, locale count and threads), the
restarted run produces *identical* results — every
subsequent dump is byte-for-byte equal to the uninterrupted run's.
This is enforced in CI:

- `restart-sod` — a Sod run is stopped after its first step (a
  pre-planted `stop` file), resumed, and every dump is required to be
  **identical** to an uninterrupted reference run;
- `restart-turbulence` — the same protocol on the OU-driven turbulence
  box with tracer particles and periodic `--restartDt` dumps: the
  hardest case, since it also requires the random forcing sequence and
  the particle trajectories to continue exactly.

Caveats worth knowing:

- a recompile (different compiler version, flags, or optimisation
  level) may change floating-point instruction scheduling — the resumed
  run is then still *correct*, just not bit-identical;
- restarting is exact because the `stop` file interrupts *between*
  steps. Runs that ended at `--tstop` clamped their last `dt` to hit
  the stop time exactly, so continuing a *finished* run to a larger
  `--tstop` is a new trajectory after that point (correct, but not the
  same dt sequence as one long run would have taken);
- restart writes rotate `restart.chaa` → `restart.bak.chaa`, so the
  two most recent checkpoints are always on disk; copy a file aside if
  you want to keep an older one.
