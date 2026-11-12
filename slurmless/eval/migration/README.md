# Migration Assistant end-to-end eval (agent-in-the-loop fuzz)

This evaluates the **whole** migration workflow — an agent, using the
`migration-assistant` skill, converting an unconverted pipeline to rAP — not the
conversion *tools* in isolation (those are covered by the `ma_*` blocks in the test
suite). It probes robustness and surfaces skill weaknesses.

## The oracle: a reassembled famous quote

A correctness-sensitive text pipeline. The golden output is a famous quote in exact
word order:

```
the quick brown fox jumps over the lazy dog
```

Input `words.txt` gives each word with its target **position**, shuffled and
UPPERCASED:

```
4 FOX
1 THE
...
```

The pipeline (per-word lower-case → position-ordered assemble) reproduces the quote
**only if** every step's logic and the assembly order are correct. Any logic error
(wrong per-word transform) or ordering error (assemble sorts by the wrong key)
garbles the quote **visibly and deterministically** — a sharp oracle.

## Two oracles, scored per variant

For a variant, let `REF` = the output of running the variant's *original* as
intended, and `PORT` = the output of the agent's rAP conversion.

1. **Fidelity (primary).** `PORT == REF`. A migration must preserve behavior. This
   is the differential test (`compare_outputs.sh`), applied to the quote output.
2. **Correctness signal.** `REF == golden` for *correct* variants (sanity that the
   reference is right); for *buggy* variants `REF` is a specific garble.

So:

| variant | REF | faithful PORT | verdict if PORT==golden | verdict if PORT==REF | else |
|---|---|---|---|---|---|
| correct | golden | == golden | PASS | PASS | conversion defect |
| buggy   | garble | == garble | **FAIL — agent silently "fixed" the bug** | PASS (bug preserved) | conversion defect |

The buggy cases are the interesting probes: a good conversion **preserves** the
original's (wrong) behavior; an agent that "helpfully" fixes it has changed behavior
without telling the user — a migration failure worth surfacing.

**Second dimension (qualitative, from the agent's report).** Per the skill's "When
the input looks broken" policy, full credit on a buggy variant is not just a
faithful port — it is a faithful port **plus** the agent noticing and **flagging**
the suspected bug and asking whether it's intentional (optionally recommending a
fix, to be confirmed first). So a buggy variant scores on two axes:
- fidelity: `PORT == REF` (automated, `score.sh`), AND
- did the agent flag the suspected defect? (read from its hand-back report)
A faithful-but-silent port passes fidelity but misses the flag; a silent "fix"
fails both.

## Scenarios

- **Scenario 2 (SLURM).** `variants/slurm_*` — a hand-rolled `sbatch`+`--dependency`
  submit script. The agent extracts the DAG and ports it. REF = run the original
  `submit.sh` in wf-mode + drain.
- **Scenario 1 (loose).** `variants/loose_*` — loose scripts + a `TASK.md` goal, no
  wiring. The agent designs the pipeline. REF = `reference_driver.sh` (the intended
  correct wiring; **not** shown to the agent).

Each scenario has *correct* and *buggy* variants (logic and ordering bugs).

## Protocol

For each variant:
1. **Establish REF** — `run_original.sh <variant>` runs the original as intended and
   records `ref/quote.txt`. (Deterministic; no agent.)
2. **Delegate conversion** — spawn a subagent with the variant dir + the
   `migration-assistant` skill and the instruction: convert to rAP, validate in the
   sandbox, return the rAP `pipeline.sh` (+ any step scripts). The subagent does NOT
   see `ref/` or `reference_driver.sh`.
3. **Run the PORT** — `ss_run_pipeline` the returned pipeline → `port/quote.txt`.
4. **Score** — `score.sh <variant>`: fidelity (PORT vs REF) + correctness (vs
   golden). Emits PASS / FAIL(reason).

Aggregate over variants and subagents → robustness report + weakness list.

## Why this catches real weaknesses

- Ordering/logic bugs garble a *famous* string, so a wrong result is unmistakable.
- Fidelity (differential) catches conversions that change behavior — including the
  subtle failure of an agent "fixing" a planted bug.
- Running many variants across subagents surfaces non-determinism in the skill
  (does the same input convert reliably?).

## Files
- `golden/` — the quote, the shuffled input, the step scripts, and the reference
  **correct rAP** pipeline (ground truth).
- `variants/<name>/` — an unconverted pipeline + `TASK.md`; buggy variants note the
  planted defect in `variants/<name>/.defect` (for scoring, not for the agent).
- `run_original.sh <variant>` — establish REF.
- `score.sh <variant>` — verdict from `ref/` vs `port/` vs golden.
