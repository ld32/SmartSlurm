# rAP annotation cheat-sheet

`runAsPipeline` converts an annotated bash script into a SmartSlurm pipeline. A step
is a marker line + its body + a close:

```
#@<stepNum>,<dependsOnStep>,<stepName>,,<inputsVar>,<sbatch options>
<command(s) — copied verbatim from the original, in the user's style>
#@end
```

Field by field:
- `stepNum`      — this step's number (1, 2, 3, …).
- `dependsOnStep`— the step number this one depends on (`0` = none). Fan-in: a
  downstream step naming an upstream step that fanned out waits for all its jobs.
- `stepName`     — short name (becomes the job name / flag).
- (empty field)  — reserved.
- `inputsVar`    — the shell var holding this step's input path(s); drives resource
  estimation and the outOfRange guard. `none` for input-less steps.
- `sbatch options` — e.g. `sbatch -p short -c 1 --mem 4G -t 50:0`.

Fan-out / loops:
- A `for` loop **around** a `#@ … #@end` block makes one job per iteration
  (per-sample). The loop variable(s) are appended to the job flag so jobs don't
  collide.
- A nested `for` inside produces per-chunk fan-out. Do **not** assign a loop var
  *inside* a block and expect it to be node-scope — assign it in the enclosing loop.

Block closing:
- Explicit: `#@end` (recommended; the whole script then uses explicit mode).
- Legacy: a blank line closes a block (fragile around loop `done` — prefer `#@end`).

Minimal example (per-sample chain + fan-in):
```sh
for s in $samples; do
    in=doc.$s.txt
    #@1,0,normalize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr 'A-Z' 'a-z' < $in > norm.$s.txt
    #@end
    #@2,1,tokenize,,in,sbatch -p short -c 1 --mem 100M -t 5:00
    tr ' ' '\n' < norm.$s.txt > tokens.$s.txt
    #@end
done
#@3,2,merge,,in,sbatch -p short -c 1 --mem 100M -t 5:00
cat tokens.*.txt > all.txt
#@end
```

Run it: `runAsPipeline --script "pipeline.sh <args>" --sbatch-options "<opts>" --tmp noTmp`
(`--mode dryrun` previews the conversion; `--lint` checks it without submitting).
