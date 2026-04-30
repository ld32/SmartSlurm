#!/usr/bin/env bash
#
# diag_stuck.sh PID
#
# Quick diagnosis of a "stuck" process:
#  - short strace to find the main blocking syscall
#  - if wait4(...) on a child, recursively inspect the child
#  - if read(...) on fd N, inspect N via lsof (/proc)
#
# Requires: strace, lsof, awk, grep, sed, head

set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <PID>" >&2
    exit 1
fi

ROOT_PID="$1"

if ! kill -0 "$ROOT_PID" 2>/dev/null; then
    echo "PID $ROOT_PID does not exist or is not accessible." >&2
    exit 1
fi

TMPDIR="${TMPDIR:-/tmp}"
TRACE_TIME=5     # seconds of strace sampling

trace_it() {
    local pid="$1"
    local indent="${2:-}"

    if ! kill -0 "$pid" 2>/dev/null; then
        echo "${indent}PID $pid no longer exists."
        return
    fi

    echo "${indent}=== Diagnosing PID $pid ==="
    ps -fp "$pid" | sed "1d" | sed "s/^/${indent}/"
    echo

    local trace_file="$TMPDIR/diag_stuck_${pid}_$$.strace"

    echo "${indent}Running strace on $pid for ${TRACE_TIME}s..."
    # -ttt: absolute timestamps, -f follow threads, -s 200 capture up to 200 chars
    # Run strace in the background; kill it after TRACE_TIME.
    ( strace -p "$pid" -f -s 200 -ttt -o "$trace_file" ) &
    local spid=$!
    sleep "$TRACE_TIME" 2>/dev/null || true
    kill "$spid" 2>/dev/null || true
    wait "$spid" 2>/dev/null || true

    if [[ ! -s "$trace_file" ]]; then
        echo "${indent}No strace output collected (permission issue or process exited?)."
        return
    fi

    echo "${indent}Sample strace (tail):"
    tail -n 10 "$trace_file" | sed "s/^/${indent}    /"
    echo

    # Find the most frequent *blocking* syscall from the tail (last ~100 lines).
    # This is a heuristic; adjust as needed.
    local top_syscall
    top_syscall="$(tail -n 100 "$trace_file" \
        | sed 's/^[0-9]\+\s\+[0-9]\+\s\+//' \
        | sed 's/^\([a-zA-Z0-9_]\+\).*/\1/' \
        | grep -E 'read|write|select|poll|ppoll|epoll_wait|wait4|futex|nanosleep|clock_nanosleep' \
        | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}' || true)"

    if [[ -z "$top_syscall" ]]; then
        echo "${indent}Could not identify a dominant blocking syscall from recent trace."
        return
    fi

    echo "${indent}Dominant blocking syscall: $top_syscall"

    case "$top_syscall" in
        wait4)
            # Look at which child(ren) it's waiting on.
            echo "${indent}Inspecting children (wait4)..."
            ps --ppid "$pid" -o pid,ppid,user,cmd | sed "1d" | sed "s/^/${indent}    /"
            local child
            child="$(ps --ppid "$pid" -o pid= | head -n1 | awk '{print $1}' || true)"
            if [[ -n "$child" ]]; then
                echo "${indent}Recursively diagnosing child PID $child:"
                echo
                trace_it "$child" "${indent}    "
            else
                echo "${indent}No children found, but wait4 is dominant (maybe zombie child or race)."
            fi
            ;;
        read)
            # Find the fd most commonly used in read(...) from the tail.
            local fd
            fd="$(grep 'read(' "$trace_file" | tail -n 50 \
                  | sed 's/.*read(\([0-9]\+\).*/\1/' \
                  | sort | uniq -c | sort -nr | head -n1 | awk '{print $2}' || true)"
            if [[ -z "$fd" ]]; then
                echo "${indent}read() is dominant, but could not determine FD."
                return
            fi
            echo "${indent}read() appears to be blocking on fd $fd"
            echo "${indent}lsof for fd $fd:"
            lsof -p "$pid" 2>/dev/null | awk -v fd="$fd" '$4 ~ ("^"fd) {print}' \
                | sed "s/^/${indent}    /" || true

            # If it's a FIFO/pipe, try to find the peer.
            local inode
            inode="$(lsof -p "$pid" 2>/dev/null | awk -v fd="$fd" '$4 ~ ("^"fd) && $5=="FIFO" {print $10}' | head -n1 || true)"
            if [[ -n "$inode" ]]; then
                echo "${indent}FD $fd is a FIFO/pipe (inode=$inode). Searching for the writer/peer..."
                # Look for other processes with same pipe inode.
                lsof -nP 2>/dev/null | awk -v inode="$inode" '$10==inode {print}' \
                    | sed "s/^/${indent}    /" || true
            fi
            ;;
        futex)
            echo "${indent}Process is blocked in futex() (likely waiting on a lock/condition within the program)."
            echo "${indent}Use a Python/Ruby/Java stack tracer (e.g. py-spy, gdb) for language-level detail."
            ;;
        select|poll|ppoll|epoll_wait)
            echo "${indent}Process is blocked in $top_syscall() (I/O wait or periodic polling)."
            ;;
        *)
            echo "${indent}Process is mostly blocked in $top_syscall(); further analysis is workload-specific."
            ;;
    esac

    echo
}

trace_it "$ROOT_PID"