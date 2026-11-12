# Task: build a SmartSlurm / runAsPipeline pipeline from these scripts

You have two step scripts and an input, but no pipeline yet:
- `lower.sh <position>` — emits `w.<position>.txt` for one word position in `words.txt`.
- `assemble.sh` — concatenates the per-word outputs in position order into `quote.txt`.
- `words.txt` — one `<position> <WORD>` per line.

Using the migration-assistant skill (Scenario 1), design and author an annotated
`runAsPipeline` pipeline that runs `lower` for every word position, then `assemble`,
and validate it in the slurmless sandbox. Leave `pipeline.sh` here. Keep the step
payloads in their original form.
