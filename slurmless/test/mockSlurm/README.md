# mockSlurm — Slurm command stand-ins

Each file here is a real executable placed **first on `$PATH`** so SmartSlurm's
unmodified calls resolve to it. Each reimplements only the slice of the real
tool's observable behavior that SmartSlurm actually depends on.

## `sbatch`
The important one. Models `sbatch [opts] --parsable -o OUT -e ERR -J NAME [-d afterok:IDS] [-H] SCRIPT`:
- prints **only the numeric job id** with `--parsable` (what `ssbatch` parses),
- **executes `SCRIPT` for real** with `SLURM_JOBID` exported, redirecting its
  output to `OUT` — so the `cleanUp.sh` EXIT trap, `memCpuMonitor.sh`, and the
  payload all run and produce genuine `.success`/`.failed` flags and `.out` logs,
- records `id`/`name`/`deps`/`hold` and the injected failure mode under
  `$MOCK_SLURM_STATE/jobs/<id>.{meta,rc,inject}` for `sacct` to replay.

Arg-parsing note: options that take a value (`-o`, `-e`, `-J`, `-d`, `-p`, `-c`,
`-t`, `--mem`, …) **skip their value token** so the value (e.g. the `.out` path)
is never mistaken for the job script. Getting this wrong silently runs the wrong
file — the first bug this harness caught.

## `sacct`
Replays recorded state in Slurm's **default whitespace-aligned** column format
(not CSV), because `cleanUp.sh` parses it positionally: `cut -d" " -f5` for
MaxRSS and `-f7` for NodeList on the **last** line, plus a substring grep of the
whole blob for the State token. Since `cleanUp.sh` reads that last line
*unquoted* (so empty fields collapse away), the mock keeps fields 1–7 non-empty
on every row. Like a real `sacct -j <id>`, it emits the allocation line plus
`<id>.batch` and `<id>.extern` step lines, so `tail -n 1` lands on `.extern`
exactly as on the cluster.

It honors `--format=<cols>` (drives which columns print, `%width` stripped),
`--noheader`/`-n` (used by `snakemakeJobStatus.sh` with `--format=State`),
`-p`/`--parsable` (pipe-delimited), and `--units`. The
`sacct --name <flag> --format jobid` lookup (ssbatch's socket-timeout fallback)
returns empty.

State = `COMPLETED` when the recorded rc is 0 and nothing was injected; otherwise
`OUT_OF_MEMORY` / `TIMEOUT` / `FAILED` / `CANCELLED` per the injected mode or a
non-zero rc. `cleanUp.sh` greps for the *partial* `OUT_OF_ME` and `CANCEL`
(because real Slurm truncates long State strings in narrow columns); the full
strings emitted here contain those substrings, so detection works either way.

> **Version pin: Slurm 25.11.7 (production).** The 25.11 release changed none of
> `sacct`'s tabular output, column names, or State enumerations, and the columns
> are fixed by `cleanUp.sh`'s explicit `--format=` string, so this format matches
> the production cluster. If you later move to a Slurm major version that changes
> `sacct` output, capture `sacct --version` and a sample row and update the
> `colval`/row emission here. This remains the single most version-sensitive file.

## `srun`
Runs a task locally. Consumes only srun's **leading** options (those before the
executable); once the executable is reached, the rest of the command line is
passed through **verbatim**. This matters because the job body is
`srun -n 1 sh -e -o pipefail -c '…'` — the `-e -o pipefail -c` belong to `sh`,
not srun, and must not be stripped (doing so was the second bug this harness caught).

`srun` is also where **failure injection** lands: when `mock sbatch` marks a job
for injection it exports `MOCK_INJECT`, and `srun` then emits a representative
scheduler message (cgroup OOM kill, `DUE TO TIME LIMIT`, or a non-zero task exit)
and exits non-zero **before** running the payload — so the job body's
`&& touch .success` never fires and `cleanUp.sh` falls through to the `sacct`
state, exactly like a real kill. This is why an injected OOM job ends up with no
`.success` flag and `jobStatus: OOM`, rather than a spurious success.

## `squeue`, `scontrol`, `sshare`, `sinfo`
Minimal. Nothing is ever pending/running (synchronous model), so `squeue` is
empty; `scontrol requeue|hold|release|update` are accepted no-ops; `sshare`/
`sinfo` are stubs. Extend them only when a code path you're testing calls them
in a way that needs a real answer.

## `scancel`
Marks the given job ids as `cancelled` in mock state.

## State directory
Everything is recorded under `$MOCK_SLURM_STATE` (set per-workspace by the
harness). Delete it to reset. It never touches a real Slurm or `~/.smartSlurm`.
