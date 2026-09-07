# Pipeline Orchestration: YAML-Driven Automated Chain

Runs Scripts 1 through 8 as one real, unattended Slurm job chain instead
of submitting and checking each one by hand. This uses **Slurm's own
native job dependency mechanism** (`--dependency=afterok:<jobid>`) -- no
external workflow-manager software (Snakemake, Nextflow, etc.) is
installed or required. `pipeline_config.yaml` just declares the step
order and script names; `run_pipeline_from_yaml.sh` reads it and submits
the whole chain.

## Quick start

```bash
cd pipeline_orchestration
bash check_environment.sh && bash run_pipeline_from_yaml.sh pipeline_config.yaml
```

That's it -- this returns immediately after submitting all 8 jobs. Slurm's
scheduler runs each one in turn, only starting step *N+1* once step *N*
has actually exited successfully.

## Preflight check (`check_environment.sh`)

A real orchestrated run on a sibling repo in this curriculum
(PRJNA274277) failed partway through a chain because a script's hardcoded
module version didn't exist on this cluster's module tree -- that only
surfaced after earlier steps had already run and consumed real job time.
`check_environment.sh` catches this whole category of failure *before*
submitting anything real:

- checks every `module load` referenced in the pipeline against what
  `module spider` reports this cluster actually has;
- compares the manual walkthrough scripts
  ([`../pipeline_raw_reads_to_counts/`](../pipeline_raw_reads_to_counts/))
  against the orchestrator's own copies ([`scripts/`](scripts/)) to catch
  drift -- since a fix has to be applied to both copies by hand, this
  check exists specifically to catch the case where it was only applied
  to one;
- checks that every referenced reference-genome download URL is
  reachable.

```bash
bash check_environment.sh
```

Exits 0 if everything checks out, 1 otherwise, so it's safe to chain
directly in front of the real run as shown above. It does not -- and
cannot -- auto-fix anything; it only surfaces "the environment doesn't
have what this script assumes it has" before you burn a job slot finding
out the hard way.

## How the chaining actually works

Each step is submitted with `sbatch --parsable`, which prints back just
the new job's ID. That ID is fed into the *next* step's submission as
`--dependency=afterok:<that_id>`. Every script's own real `#SBATCH`
resource directives (`--mem`, `--cpus-per-task`, `--time`, etc.) are left
completely untouched -- the orchestrator only adds the dependency flag on
top.

`8_QC_Report_Generator.sh` is documented in its own header to run locally
via plain `bash`, not `sbatch` -- for manual use, keep doing that. For the
automated chain specifically, it's wrapped in a minimal 1-CPU/2GB/10-minute
job (`mode: sbatch_wrap` in the YAML) purely so Slurm has something to
attach the final dependency to. This doesn't change what the script does,
only how the chain triggers it.

## Isolation from manual runs: separate script copies, not shared ones

This orchestrator does **not** touch the student-facing walkthrough at
all:

- [`../pipeline_raw_reads_to_counts/`](../pipeline_raw_reads_to_counts/) -- the student-facing walkthrough. Byte-identical to before this orchestrator existed. A student runs each script by hand exactly as the main [`README.md`](../README.md) describes, and nothing here changes that experience.
- [`scripts/`](scripts/) -- **copies** of the same 8 scripts, used *only* by this orchestrator. Each copy's `WORKDIR` gets an always-distinct suffix:

```bash
WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep_${RUN_TAG:-orchestrated}"
```

Because the fallback (`orchestrated`) applies even if `RUN_TAG` is
somehow left unset, an orchestrated run's working directory can never
collide with a manual run's directory -- there's no shared "untagged"
path.

Why two copies instead of parameterizing the manual scripts in place?
So a student who opens a script in `pipeline_raw_reads_to_counts/` to do
the walkthrough by hand never has to understand or even see the
`RUN_TAG` mechanism -- the teaching script stays exactly as simple as it
was designed to be. The trade-off, made deliberately: if you fix a bug in
one copy, you must copy the fix to the other -- `check_environment.sh`'s
drift check above exists specifically to catch when that's forgotten.

Any shared reference-genome/index directory this pipeline builds is
**not** suffixed in either copy -- it's expensive to build and every
script that builds it already checks for its output before rebuilding,
so sharing it across manual and orchestrated runs alike is safe and
avoids redundant work.

Want a second, independent orchestrated run without deleting the first
one's output (e.g. comparing two configurations)? Give it a different
tag:

```bash
bash run_pipeline_from_yaml.sh pipeline_config.yaml --run-tag demo2
```

## Monitoring progress

```bash
squeue -u $(whoami)
```

Shows every step currently queued or running. A step waiting on its
dependency shows `(Dependency)` as its reason.

## If a step fails

Slurm's `afterok` dependency means: if step *N* exits non-zero, every
step after it is left in a `DependencyNeverSatisfied` state and **never
runs** -- it does not retry, skip, or silently continue with wrong data.
Find out exactly which step failed and why:

```bash
sacct -u $(whoami) --format=JobID,JobName%20,State,ExitCode,Start,End
```

Fix whatever caused that step to fail in `scripts/` (the orchestrator's
own copy), then **resume from that exact step** instead of re-running the
whole chain from scratch:

```bash
bash run_pipeline_from_yaml.sh pipeline_config.yaml --start-from 3
```

## Editing the pipeline

`pipeline_config.yaml`'s `steps:` list is deliberately written as
one-line YAML flow mappings (`- {id: N, script: ..., mode: ...}`) so it
can be parsed reliably with plain `grep`/`sed` on any HPC login node --
no `PyYAML` or `yq` installation required. It's still fully valid YAML;
verify it any time with:

```bash
python3 -c "import yaml; print(yaml.safe_load(open('pipeline_config.yaml')))"
```

To add, remove, or reorder a step, edit that list -- `mode: sbatch` for a
script with its own real `#SBATCH` header, `mode: sbatch_wrap` for a
script meant for local/manual use that you still want inside the
automated chain. If you change a script's real logic, remember to make
the same change in both `../pipeline_raw_reads_to_counts/` (student copy)
and `scripts/` (orchestrator copy) -- see "Isolation from manual runs"
above for why they're kept separate.
