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

## User-workflow mode (`SS_WF_MODE`)
Everything above is the **default** mode (jobs run synchronously in submission
order — how SmartSlurm's own tools are tested). When `ss_wf_activate` sets
`SS_WF_MODE=1`, three commands switch to a **register-then-drain scheduler** so an
**unmodified** hand-rolled SLURM workflow runs faithfully (this is what the
Migration Assistant builds on). Default mode is unaffected.

- **`sbatch` — register, don't run.** Instead of executing the script, it records a
  job node under `$MOCK_SLURM_STATE/wf/<id>.{name,deps,out,err,kill,array,run}`
  (command, parsed `--dependency`, `--array`, `--wait`, `--kill-on-invalid-dep`,
  resources), prints `Submitted batch job <id>`, and **returns immediately** — so
  the user's `jid=$(sbatch …)` capture + later `--dependency=afterok:$jid` build a
  real DAG. The job id comes from an **atomic** `.jobseq` counter (flock, with an
  mkdir-spinlock fallback) so concurrent submits never collide. `--wait` drains that
  job and its deps inline and exits with the job's rc.
- **`wfdrain` — the DAG executor.** Not a Slurm command; invoked by `ss_wf_drain`
  after the submit script returns (the natural "submit-all" barrier). It evaluates
  each job's dependency expression to READY / PENDING / INVALID and runs READY jobs
  to a fixpoint, honoring the **full grammar** (after / afterany / afterok /
  afternotok / afterburstbuffer / aftercorr / singleton; `,`=AND / `?`=OR; `+time`),
  arrays (per-task `SLURM_ARRAY_TASK_ID`, element-wise `aftercorr`), fault injection
  (`MOCK_FAIL_*` → rc≠0 + `.failmode` + `FAILED(mode)`), and tiered stubbing
  (`SS_WF_STUB` / `SS_WF_STUB_MISSING`). A job whose dependency can never be
  satisfied is killed or left held per `--kill-on-invalid-dep` / `SS_WF_KILL_INVALID_DEP`.
- **`scancel` / `squeue` — live state.** In wf-mode `scancel <id>` settles a pending
  job to cancelled (and cascades to dependents), and `squeue` reports the active /
  held jobs (`-h`, `-j`) instead of being empty. Default-mode behavior (above) is
  preserved when `SS_WF_MODE` is unset.

See `test/README.md` → *User-workflow mode (Migration Assistant)* for the harness
API (`ss_wf_activate` / `ss_wf_drain` / `ss_run_workflow`) and the `wf_*` test blocks.

## State directory
Everything is recorded under `$MOCK_SLURM_STATE` (set per-workspace by the
harness). Delete it to reset. It never touches a real Slurm or `~/.smartSlurm`.
In user-workflow mode the registered DAG lives under `$MOCK_SLURM_STATE/wf/`.
