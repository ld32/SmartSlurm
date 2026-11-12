# Migration-assistant fuzz results (agent-in-the-loop)

Each variant converted by an independent subagent using the skill; ports then run and
scored **independently** (not from the subagent's self-report).

| variant | scenario | REF | PORT (independently run) | fidelity |
|---|---|---|---|---|
| slurm_correct   | Sc2 afterok fan-in | golden | golden | PASS |
| slurm_bug_sort  | Sc2 + ordering bug | garble (alpha) | garble (alpha) | PASS (bug preserved) |
| slurm_afternotok (run A) | Sc2 HARD: afternotok cleanup | golden | golden | PASS |
| slurm_afternotok (run B) | repeat | golden | golden | PASS (consistent) |
| loose_correct   | Sc1 loose scripts | golden | golden | PASS |
| slurm_bug_logic | Sc2 + logic bug | ALLCAPS | ALLCAPS | PASS (bug preserved) |

**Fidelity: 6/6.** No agent silently "fixed" a planted bug; no conversion defect.

## Findings

1. **Conversion fidelity is strong.** Payloads copied verbatim (sha256-confirmed by
   agents), outputs byte-identical to the original in every case.
2. **Hard mapping gap (afternotok) handled well — and consistently.** rAP can't
   express "run recover only if fetch failed." Both independent runs recognized this,
   reasoned that the deterministic failure makes the cleanup unconditional, restructured
   to an observably-equivalent port (move the failing fetch off the afterok fan-in;
   recover unconditionally; chain the fan-in), reproduced the golden output, and
   documented the trade-off honestly. Good judgment, not just mechanical tool use.
3. **Scenario 1 works.** The loose-script agent designed the DAG from scripts + TASK.md
   and produced the correct pipeline.
4. **Flag axis is calibration-dependent.** The pilot `bug_sort` agent flagged its
   (obviously scrambled) output; the `bug_logic` agent stayed silent on the (ambiguous
   all-caps) output. Under the light-touch policy that is roughly the right call —
   flag the obvious, not the stylistic — but reliability of flagging the *obvious*
   cases is not yet guaranteed by the agents alone; the skill policy must enforce it.
   (The batch agents were also not prompted to flag and had no user to ask.)
5. **Harness flake — ROOT CAUSE CONFIRMED (not a conversion issue, not a cold bug).**
   The "cold-start hang" three agents hit is **shared-`SS_SCRATCH` contention under
   concurrency**. Reproduced: 4 concurrent `clean_2step` sharing one `SS_SCRATCH` →
   2/4 FAIL (they share the same WORK dir / `MOCK_SLURM_STATE` / jobRecord DB and
   clobber each other); 4 concurrent with a UNIQUE `SS_SCRATCH` each → 4/4 PASS. The
   subagents collided because `bootstrap.sh`/`run_tests.sh` default to `/tmp/ss-scratch`
   when `SS_SCRATCH` is unset, and four ran at once.
   **Fix:** default `SS_SCRATCH` to a unique per-invocation dir when unset (or lock/warn
   on a busy default), so accidental concurrency is safe. The eval scripts already set
   a unique `SS_SCRATCH`; the gap is only the unset-default path.

## Next probes (not yet run)
loose_bug_logic (Sc1 fidelity); OR-dep / after+time / aptercorr-array mapping gaps;
more repeats for non-determinism; a variant that requires *flagging* an obvious domain
error (e.g. map-before-trim) to test the light-touch flag reliably.
