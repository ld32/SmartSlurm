[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)

# FAQ and troubleshooting

## ssbatch FAQ

<details>
<summary><b>Do I have to wait for the first jobs to finish before later jobs get estimated resources?</b></summary>

For **standalone ssbatch**: there is no waiting — each job is submitted immediately. Estimation simply kicks in once ≥3 records exist. For **runAsPipeline**: yes, by design it runs the first `firstBatchCount` (5) independent jobs and holds the rest until estimates are available.
</details>

<details>
<summary><b>Are -P, -I, and -F optional?</b></summary>

Yes, all three.
- `-P` omitted → program name is taken from the wrapped command or the script name.
- `-I` omitted → estimation uses the program's past-usage distribution (90th percentile) instead of a size-based fit.
- `-F` omitted → the unique flag becomes `program`+`input`, or `program`+a random suffix.
</details>

<details>
<summary><b>Can -I take a size directly instead of a file?</b></summary>

Yes: `-I jobSize:12` — `12` is the input size (any integer). Useful when the meaningful "size" isn't a single file.
</details>

<details>
<summary><b>Can I pass -c and other sbatch options?</b></summary>

Yes. Any standard sbatch option works and is passed through.
</details>

<details>
<summary><b>How does ssbatch know a job already ran?</b></summary>

On success it creates `<flag>.success` in `smartSlurmLogDir`. That file's presence is how reruns are detected.
</details>

## runAsPipeline FAQ

<details>
<summary><b>Do later jobs wait for the first jobs before getting estimated resources?</b></summary>

If a step's jobs are independent, `runAsPipeline` submits them all at once, lets the first `firstBatchCount` (5) run, and holds the rest. Once early jobs finish and produce ≥3 records, the held jobs are released with estimated resources.
</details>

<details>
<summary><b>Should I use `#@end` or the blank-line close?</b></summary>

Prefer **`#@end`** for anything beyond a trivial one-liner: the body is written verbatim to a per-job script, so `awk`, quotes, heredocs, multi-line bodies, and command substitution all work. The legacy blank-line close flattens the body into a single `--wrap` string, which mangles quoted whitespace and multi-line bodies. A pipeline is all-one-mode: if it uses `#@end` anywhere, every block must close with `#@end`.
</details>

<details>
<summary><b>Can inputs be given as a size instead of a file?</b></summary>

Yes. Assign a variable and reference it in the marker:
```bash
someVar=jobSize:12   # 12 is the input size (any integer)
#@1,0,runshard,,someVar,sbatch -p short -c 4 -t 0-12:00 --mem 8G
```
</details>

<details>
<summary><b>Multiple inputs?</b></summary>

Yes — dot-join them in the marker (`#@2,1,find,,input1.input2,...`) or use a shell variable holding a space-separated list.
</details>

<details>
<summary><b>Fewer or no emails?</b></summary>

Add `noSuccEmail` (failures only) or `noEmail` (none) with `--email`.
</details>

<details>
<summary><b>Can I drop the command-line sbatch options?</b></summary>

Yes, **if every step sets its own** `sbatchOptions`. Then omit `--sbatch-options`.
</details>

<details>
<summary><b>Does runAsPipeline run my script's commands in their original order?</b></summary>

Plain commands run top to bottom during submission. Commands under a `#@` marker do **not** run then — they are submitted as Slurm jobs and run later on compute nodes (respecting dependencies). See [the two phases](core-concepts.md#how-a-pipeline-runs-convert-and-submit-two-phases).
</details>

<details>
<summary><b>How does runAsPipeline detect an already-completed step?</b></summary>

By the `<flag>.success` file in `smartSlurmLogDir` (set in `config.txt`).
</details>

---


---

[&larr; SmartSlurm README](../README.md) &middot; [Docs index](README.md)
