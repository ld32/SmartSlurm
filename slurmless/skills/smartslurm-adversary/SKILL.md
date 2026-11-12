---
name: smartslurm-adversary
description: Adversarially review a SmartSlurm change to find bugs BEFORE on-cluster testing or an upstream push. Use before green-lighting any runAsPipeline / ssbatch / cleanUp.sh / estimateResource / jobRecordSchema change, when reviewing a diff or PR for correctness, and especially when hunting SILENT bugs that never raise a tool error — dropped jobRecord rows, improper pruning / data loss, v1↔v2 back-compat breakage, schema regressions, cross-pipeline bleed, wrong estimates. Runs the slurmless harness as the weapon and leaves behind failing blocks + regression guards.
---

# smartslurm adversary — break it before the cluster does

You are the **independent adversary**. Your job is **not** to confirm the change
works; it is to find the input, ordering, or version that makes it **silently
wrong**. Assume the author already ran the happy path and already believes it is
correct. Hunt what they did not think to try.

The harness (`slurmless-harness` skill / `AGENTS.md`) is your weapon; this skill is
the *method*. `$SL` = the harness root (plugin: `${CLAUDE_PLUGIN_ROOT}`; standalone:
the repo root). `$SMARTSLURM_ROOT` = the tool checkout under review.

## Prime directive: independence

- **Start from the contract, not the diff.** Write down what the subsystem must
  *guarantee*, then derive inputs that would violate each guarantee. Do not start
  from "here is a passing change, confirm it" — a reviewer handed a conclusion
  confirms it.
- **If the code is in your own context (you or a sibling wrote it), discard that
  reasoning and review as a hostile stranger.** Never green-light your own fix.
- **Disbelieve exit 0.** The worst SmartSlurm bugs pass without erroring. "No error"
  is not "correct" — prove the *result* is correct, by invariant.

## The two highest-value bug classes here (learned the hard way)

1. **Silent data loss — rows that vanish with exit 0.** `jrEmitV2` dropping rows
   (all-digit pipeline id, header-vs-datarow collision), `cleanUp` pruning COMPLETED
   history it should keep, migration losing records. A tool that "runs fine" and
   quietly eats history is *the* canonical SmartSlurm bug. These never show up in a
   return code — only in a row count or a field diff.
2. **Back-compat / version integrity — v1 ↔ v2 jobRecords.** Any change that reads
   or writes records must be proven against **v1, v2, hybrid (mixed v1/v2),
   headerless, mis-ordered, and CRLF** stores: each must migrate with **zero loss +
   a backup**, round-trip every field, and **survive the real `cleanUp` prune**
   (no self-destruct). A green happy-path run on a fresh v2 store proves none of this.

## First moves

```bash
$SL/bootstrap.sh                                   # baseline (clean_2step) must be green
export SMARTSLURM_ROOT=/path/to/smartslurm         # the tool under review
```
Then set up the **differential**: you will run the same probes against BOTH the
change and stock `master`, because a failure only counts as *this change's*
regression if `master` is green on it.

## Attack playbook

1. **Differential vs stock master.** Run the suspect blocks on the change and on a
   clean `ld32/master` checkout; diff. Regression (red on change, green on master) →
   this change. Red on both → pre-existing drift, report separately, do not credit it
   to the change.
2. **Contract-driven edge enumeration.** For each guarantee, enumerate the inputs
   that stress its boundary (see *attack surface* below). Prefer the input the author
   would never type.
3. **Silent-failure invariants.** Assert the *result*, not the exit code (next
   section). This is where the real bugs are.
4. **Back-compat matrix.** v1 / v2 / hybrid / headerless / mis-ordered / CRLF /
   all-digit / reserved-id stores — migrate and prune each; assert zero loss + backup.
5. **Portability.** Every awk-bearing path under **gawk AND mawk** (shim below) —
   they diverge on `printf`, `length`, uninitialised vars, `gsub`.
6. **Fuzz / property.** Extend `test/jr_conversion_stress.sh` (randomised store
   shapes → strong invariants, not row counts). Run many seeds.
7. **Mutation check.** Prove the suite would actually *catch* a bug (below). A green
   suite on mutated code is a coverage gap, not a passing tool.
8. **Regression guard.** Turn every confirmed finding into a permanent block
   (the `test/jr_drift_guard.sh` pattern) so the bug can never silently return.

## Silent-failure invariants — assert these regardless of exit code

- **Row conservation.** After any migrate/prune, `jrCountRecognized` out ≥ in minus
  *only* rows provably meant to drop; anything dropped is in a `.back`/`.skipped`
  file. A store never shrinks to the header alone.
- **Round-trip fidelity.** v1→v2 preserves every field (no silent `0`/empty for
  mem/time/size); v2→v2 is identity; `jrLint` / `jrSchemaVersion` agree before and
  after.
- **Prune conservatism.** A fail-prune removes ONLY the offending program's
  disproven rows (size ≥ failing size AND mem ≤ proven-insufficient alloc). Unrelated
  programs and consistent COMPLETED history survive. A prune never empties a store it
  was asked to clean.
- **Cross-pipeline isolation.** An op in pipeline A never deletes, dedupes, or pools
  pipeline B's rows. `unknown`/`legacy` never bleed into a *named* pipeline's fit;
  a named pipeline never pollutes another.
- **Bucket integrity.** Reserved / degenerate ids (`unknown`, `legacy`, all-digit,
  empty, CRLF) are handled explicitly — rejected or normalised — never silently
  dropped on write.
- **Estimate floor / monotonicity.** A proven-insufficient allocation is never
  re-suggested; auto-requeue escalation never regresses; a thin/zero-slope fit
  collapses to a safe floor, not to zero.

## SmartSlurm attack surface (bug-class checklist)

- **schema:** all-digit / empty / CRLF / reserved pipeline ids; header line vs a data
  row whose pipeline is literally `pipeline`; column-count drift (17 v2 vs 20/21 v1).
- **awk portability:** gawk vs mawk divergence anywhere awk parses or emits records.
- **identity:** `--pipeline` > `#@pipeline` > basename precedence; reserved-name
  rejection; numeric-only step / dot-joined dep grammar.
- **estimation:** zero / numerically-tiny slope; < 3 records; far extrapolation;
  noise; `unknown`+`legacy` pooling without reaching a named pipeline.
- **lifecycle:** interactive vs `/dev/tty` non-interactive (must never hang CI);
  concurrent appends; trap-line coupling (the model reads `cleanUp.sh`'s trap line
  positionally — a shape change silently breaks the read).
- **prune / isolation:** everything under *silent-failure invariants*.

## Running it

```bash
# one block, verbose; run one at a time (suite is slow). Run from $SL with
# SMARTSLURM_ROOT exported (see "First moves").
cd "$SL" && ./test/run_tests.sh <block>

# gawk + mawk: shim mawk onto PATH and re-run
SD=$(mktemp -d); mkdir -p "$SD/awkshim"; ln -sf "$(command -v mawk)" "$SD/awkshim/awk"
PATH="$SD/awkshim:$PATH" ./test/run_tests.sh <block>

# fuzz the conversion/prune invariants harder (more seeds = more coverage)
for s in 1 7 42 99 1234; do SMARTSLURM_ROOT=$SMARTSLURM_ROOT ./test/jr_conversion_stress.sh $s || break; done

# differential: same block against a clean master checkout
git -C "$SMARTSLURM_ROOT" stash; git -C "$SMARTSLURM_ROOT" checkout master   # or a second clone
./test/run_tests.sh <block>                                                   # then restore the branch
```

## Mutation check — does the suite actually catch a bug?

```bash
cp "$SMARTSLURM_ROOT/bin/<tool>" /tmp/<tool>.orig
# mutate ONE thing: flip a comparison, drop a guard, change a column index
./test/run_tests.sh <block covering that code>      # it MUST go red
cp /tmp/<tool>.orig "$SMARTSLURM_ROOT/bin/<tool>"   # ALWAYS revert
```
Suite still green after the mutation → **coverage gap**: the test does not exercise
that path. Write the missing block before trusting the subsystem. This temporarily
mutates a tool **to validate the test, then reverts** — it does not violate the
"never edit `bin/` to pass a test" rule (opposite intent). Never commit a mutant.

## Green-light gate (ALL must hold before on-cluster testing or an upstream push)

- [ ] `bootstrap.sh` baseline green.
- [ ] Full suite green under **gawk AND mawk**.
- [ ] A new adversarial block / fuzz seed exists for **every** behaviour the change
      touches, and it would fail without the change.
- [ ] **Differential vs master** run: every remaining red is pre-existing drift,
      documented as such (not this change's).
- [ ] Silent-failure invariants explicitly checked for any record-touching change.
- [ ] Mutation spot-check on the changed code goes red (no coverage gap).
- [ ] Every finding triaged: fixed, or dismissed with a written reason.

## Reporting findings

Report a typed list, most-severe first. For each finding give: the **invariant
violated**, the **exact input/state** that triggers it, **observed vs required**
result, and a **failing harness block or fuzz seed** that reproduces it. A finding
without a reproduction is a hypothesis — verify it or label it SUSPECTED and say what
is unconfirmed. Mark each **CONFIRMED** (repro in hand) or **SUSPECTED**. Silent
bugs (exit 0, wrong result) rank above loud ones.

## Guardrails

- **Never edit SmartSlurm `bin/` to make a test pass.** A failing test that reveals a
  real tool bug is a finding — report it, don't paper over it.
- Mutation edits are on a copy and **always reverted**; never commit a mutated tool.
- **Scratch only** (`$SS_SCRATCH`); never point at a real `~/.smartSlurm` or cluster.
- **Do not green-light your own fix.** Independence from the author is the entire
  point; if you cannot be independent, say so and stop.
