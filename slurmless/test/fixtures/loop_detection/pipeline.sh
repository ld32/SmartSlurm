#!/bin/sh
# Regression anchor for runAsPipeline's TEXTUAL loop detection.
#
# rAP decides loop nesting by scanning the driver as text: a line matching the
# glob `for*` pushes a loop-variable suffix onto the flag, a line matching `done*`
# pops it. It is a prefix match, not a bash parse, so two common constructs
# silently corrupt every downstream job flag:
#   1. a one-line `for … ; done` pushes but never pops (the `done` is mid-line,
#      not at line start);
#   2. a command whose first word merely STARTS WITH "for" (format, foremost,
#      foreach, …) is misread as a loop — its 2nd token becomes a bogus loop var.
#
# This fixture triggers both. When rAP's detection is hardened, steps 2 and 3
# should carry CLEAN flags (2.1.second, 3.2.third) — flip the asserts then.

#@1,0,first,,input,sbatch -p short -c 1 --mem 100M -t 5:0
echo first

# trap 1: one-line for...done (pushes .$x, never pops)
for x in 1 2 3; do echo $x; done

#@2,1,second,,input,sbatch -p short -c 1 --mem 100M -t 5:0
echo second

# trap 2: a submit-scope command whose first word starts with "for"
format_report --nope

#@3,2,third,,input,sbatch -p short -c 1 --mem 100M -t 5:0
echo third
