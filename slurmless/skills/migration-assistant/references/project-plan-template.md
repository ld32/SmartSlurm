# Project plan template (Scenario 1 — loose scripts)

Fill this in **with** the user (checkpoint 1) before authoring the rAP script. Keep
it short; it exists to pin the DAG and the per-step contract.

## Pipeline: <name>

**Goal:** <one line — what the pipeline produces from what>

**Samples / fan-out axis:** <how work is partitioned — e.g. one row per sample in
`manifest.txt>`; or "single run, no fan-out">. The fan-out list must be knowable at
submit time (from a manifest/dir listing), not from an upstream job's output.

**Environment:** <conda env.yml / container image / module list so the real tools
resolve in the sandbox>

## Steps

For each step:

| # | name | inputs | outputs | depends on | resources (mem/time/cores) | tool(s) |
|---|------|--------|---------|------------|----------------------------|---------|
| 1 | …    | …      | …       | (none)     | …                          | …       |
| 2 | …    | …      | …       | step 1     | …                          | …       |

## Acceptance checks (sanity, not just exit-0)

- Which output files must exist and be non-empty?
- Expected formats / simple invariants (row counts, headers, value ranges)?
- Any golden/reference output to diff against? (optional — turns Scenario 1 into a
  Scenario-2-style comparison)

## Test inputs

- Real inputs to downsample, or synthetic tiny inputs to generate (per step's
  declared input format).
