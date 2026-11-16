# Test fixtures — harness infrastructure (do not hand-edit)

These directories are **slurmless test fixtures**, not demos. Each is a tiny
pipeline that a `run_tests.sh` block drives through mock Slurm to assert one
specific rAP / ssbatch behavior. The automated regression suite depends on them
byte-for-byte: the asserts in `../run_tests.sh` are written against the exact
step IDs, loop structure, filenames, and output of the `pipeline.sh` here.

**Do not edit a fixture to "try something" or to learn the syntax.** Changing a
`pipeline.sh` (or its `manifest.txt`) will silently change what the matching
test block proves, or turn it red. Treat this tree as code owned by the harness.

## Want a pipeline to run, copy, or learn from?

Use the **user-facing demos**, which live separately and are meant to be edited
and copied freely:

    SmartSlurm/example-scripts/        # demo_*.sh — runnable demos + README

Those show current `#@…#@end` syntax on realistic little pipelines; this tree is
the regression net underneath SmartSlurm itself.

## Legitimately changing a fixture

Only when you are intentionally changing harness coverage:

1. Edit the `pipeline.sh` / `manifest.txt` here, **and**
2. update the matching block in `../run_tests.sh` in the same change so its
   asserts still describe what the fixture does, **and**
3. re-run that block (`./test/run_tests.sh <name>`) plus the mawk-portability
   check before committing.

See `../README.md` (“Adding a fixture”) for the full block-authoring pattern.
