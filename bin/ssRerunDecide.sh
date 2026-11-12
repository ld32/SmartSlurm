#!/bin/bash
# ssRerunDecide.sh -- decide the pipeline re-run SCOPE once, UP FRONT, before the
# Stage-2 progress table is drawn.
#
# This logic used to live inside the per-step `ssbatch` call (its pipeline-reissue
# branch), where it fired lazily on the FIRST already-completed step. That had three
# bad effects: it prompted even during a `--mode dryrun`; it interrupted the first
# progress row mid-render (orphaning that sample's glyph onto its own line without a
# label); and it wrote scope markers during a dry run. runAsPipeline now calls this
# once in the generated driver's preamble, before any row is printed. `ssbatch` is
# UNCHANGED -- it still just reads the skip*/reRun* marker files this writes, so by
# the time the step loop runs the scope is already decided and ssbatch never prompts.
#
#   Env in : smartSlurmLogDir  -- dir holding the *.success flags + allJobs.txt
#            smartSlurmRerun    -- optional non-interactive override: all|failed|fromStep:N
#   Arg 1  : "dryrun" on a dry run -> never prompt; default to 'failed'
#   Exit   : 0 proceed (markers written, or nothing to decide); 2 cancel (abort the run)
#
# Markers written (exactly the names ssbatch reads at its lines ~442/445):
#   reRunAllSuccessJobs.txt           all successful steps re-run
#   skipAllSuccessJobs.txt            all successful steps kept
#   reRunAllSuccessJobs<N>.txt        step N (and up) re-run      } from-step
#   skipAllSuccessJobs<N>.txt         step N kept                 }

logDir="${smartSlurmLogDir:-smartSlurmLog}"

# Nothing completed in a prior run -> no re-run decision to make.
ls "$logDir"/*.success >/dev/null 2>&1 || exit 0

# Completed step numbers (the steps you can re-run from).
reachable=$(ls "$logDir"/*.success 2>/dev/null | sed 's#.*/##' | cut -d. -f1 | sort -un | tr '\n' ' ')
reachable=$(echo $reachable)   # trim

if [ -n "${smartSlurmRerun:-}" ]; then
    x="$smartSlurmRerun"                       # scripts / CI / tests
elif [ "${1:-}" = dryrun ]; then
    x=failed                                   # a dry run previews the default, never prompts
else
    {
        echo
        echo "!!! Some jobs already completed. What would you like to re-run?"
        echo
        echo "all      re-run ALL successful jobs (redo completed work)."
        echo
        echo "failed   re-run only failed/incomplete jobs, keep successful ones. [default]"
        echo
        echo "step     re-run from a step onward (pick from completed steps: $reachable)."
        echo
        echo "cancel   abort without submitting."
        echo
        printf 'Choice [failed]: '
    } >&2
    read x </dev/tty || x=""
    echo "You typed: \"$x\"" >&2
fi
x=${x:-failed}

case "$x" in
    all|yall)
        touch "$logDir/reRunAllSuccessJobs.txt" ;;
    failed|nall)
        touch "$logDir/skipAllSuccessJobs.txt" ;;
    cancel|q|quit)
        echo "Cancelled -- no jobs submitted." >&2
        exit 2 ;;
    step|ystep|from|fromStep|fromStep:*|from:*)
        fromStep=""
        case "$x" in fromStep:*|from:*) fromStep=${x#*:};; esac
        if [ -z "$fromStep" ]; then
            echo "Completed steps you can re-run from:" >&2
            for s in $reachable; do
                nm=$(awk -v s=$s 'NR>1 && $3 ~ ("^"s"\\.") {split($3,f,"."); print f[3]; exit}' "$logDir/allJobs.txt" 2>/dev/null)
                echo "    $s  $nm" >&2
            done
            printf 'Re-run from which step? ' >&2
            read fromStep </dev/tty || fromStep=""
            echo "You typed: \"$fromStep\"" >&2
        fi
        if [[ " $reachable " != *" $fromStep "* ]]; then
            echo "Step '$fromStep' is not in the completed set ($reachable); nothing to re-run there." >&2
            # invalid target -> keep all successful (safe default); do not abort the run
            touch "$logDir/skipAllSuccessJobs.txt"
            exit 0
        fi
        # keep steps below the chosen one, re-run it and everything above (downstream
        # then propagates through the live dependencies).
        for s in $reachable; do
            if [ "$s" -ge "$fromStep" ]; then touch "$logDir/reRunAllSuccessJobs$s.txt"
            else touch "$logDir/skipAllSuccessJobs$s.txt"; fi
        done ;;
    *)
        echo "Unrecognized choice; keeping successful jobs (re-running failures only)." >&2
        touch "$logDir/skipAllSuccessJobs.txt" ;;
esac
exit 0
