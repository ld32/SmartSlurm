# Identity contract (input to `tools/compare_outputs.sh`)

A differential test asks: does the rAP port produce the *same* results as the
original? "Same" is rarely byte-identical — runs differ in timestamps, job ids,
absolute paths, and sometimes line ordering. The identity contract is where you
(agent) state what actually counts, and the user approves it (checkpoint 2) before
the comparison is trusted.

```
# which outputs are authoritative (globs ok; repeatable)
FILE grandtotal.txt
FILE report.txt
FILE wc.*.txt

# normalize volatile fields before diffing (ERE sed; repeatable)
NORMALIZE s/[0-9]{4}-[0-9]{2}-[0-9]{2}/<DATE>/g
NORMALIZE s/job [0-9]+/job <ID>/g
NORMALIZE s#/tmp/[^ ]+#<PATH>#g

# compare these order-insensitively (repeatable)
SORT manifest.txt
```

Guidance:
- List only the files that represent the pipeline's *result* — not logs, not
  `.success` markers, not scheduler scratch. A file present in only one run is a
  mismatch.
- Prefer the narrowest normalization that removes a *legitimately* volatile field.
  Over-normalizing can mask a real regression — if you must strip a data field to
  get a match, that is a finding to raise, not a rule to add.
- `compare_outputs.sh <origDir> <rapDir> contract` prints MATCH/DIFF per file and
  exits non-zero on any mismatch.
