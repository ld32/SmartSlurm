#!/usr/bin/env python3

# Terminal-based SLURM job record reviewer — no X11, browser, or extra packages.
# Usage: reviewJobRecords.py [jobRecord.txt]
#
# jobRecords v2: columns are addressed by NAME (csv.DictReader), and records are
# grouped by (pipeline, step, ref) -- the same key the estimator fits on -- so two
# pipelines that share a program/step stay separate. A legacy v1 file is refused
# with a pointer to jrRepair.sh (this tool is offline; a hard prompt is fine).

import csv
import time
import shutil
import os
from pathlib import Path
import glob
import sys
import tempfile

# ── File path resolution ───────────────────────────────────────────────────────

if len(sys.argv) == 2:
    job_file = Path(sys.argv[1])
    smartSlurmJobRecordDir = str(job_file.parent)
else:
    script_dir = Path(__file__).resolve().parent
    home_config = Path.home() / ".smartSlurm/config/config.txt"
    local_config = script_dir.parent / "config/config.txt"

    if home_config.is_file():
        config_path = home_config
    elif local_config.is_file():
        config_path = local_config
    else:
        raise FileNotFoundError("Config file not found: config.txt")

    smartSlurmJobRecordDir = None
    with open(config_path, "r") as config_file:
        for line in config_file:
            if line.startswith("export smartSlurmJobRecordDir="):
                smartSlurmJobRecordDir = line.strip().split("=", 1)[1].replace(
                    "$HOME", str(Path.home())
                )
                break

    if not smartSlurmJobRecordDir:
        raise ValueError("smartSlurmJobRecordDir not found in config.")

    job_file = Path(smartSlurmJobRecordDir) / "jobRecord.txt"

if not job_file.is_file():
    raise FileNotFoundError(f"Job record file not found: {job_file}")

print(f"Job file : {job_file}")

tmp_dir = tempfile.gettempdir()
tmpFile = os.path.join(tmp_dir, "jobRecord_tmp.txt")
shutil.copy(job_file, tmpFile)
print(f"Working copy: {tmpFile}")

# ── Data I/O (v2, name-addressed) ───────────────────────────────────────────────

def read_job_records(ignore_patterns=None):
    """Return (fieldnames, rows) where rows are dicts keyed by column name.
    Refuses a legacy v1 file (header must start with 'pipeline')."""
    if ignore_patterns is None:
        ignore_patterns = ["missingInputFile"]
    with open(tmpFile, newline="") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or reader.fieldnames[0] != "pipeline":
            sys.exit(
                "This jobRecord is not in the v2 schema (its header does not start "
                "with 'pipeline').\nRun  jrRepair.sh  once to convert it, then re-run "
                "this tool."
            )
        rows = []
        for row in reader:
            isz = (row.get("inputSize") or "")
            if any(p.lower() in isz.lower() for p in ignore_patterns):
                continue
            rows.append(row)
    return reader.fieldnames, rows


def group_key(row):
    """The fit key: (pipeline, step, ref)."""
    return (row.get("pipeline", ""), row.get("step", ""), row.get("ref", ""))


def group_label(key):
    return "  /  ".join(k if k != "" else "-" for k in key)


def unique_groups(rows):
    return sorted({group_key(r) for r in rows})


def write_tmp(headers, rows):
    with open(tmpFile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


def save_to_file(headers, rows):
    try:
        mtime = time.strftime("%Y%m%d%H%M%S", time.localtime(job_file.stat().st_mtime))
        backup = Path(smartSlurmJobRecordDir) / f"jobRecord_backup_{mtime}.txt"
        shutil.copy(job_file, backup)
        print(f"Backup : {backup}")
    except Exception as e:
        print(f"Warning — could not create backup: {e}")
    shutil.copy(tmpFile, job_file)
    print(f"Saved  : {job_file}")


def delete_stats_files(pipeline, step):
    """Drop the cached fit stats for this pipeline+step (all ref groups) so the
    estimator refits from the trimmed history. Stats are keyed pipeline.step.refGroup.*
    (extraMem is step-keyed and pooled across pipelines)."""
    stats_dir = os.path.join(smartSlurmJobRecordDir, "stats")
    for pat in [f"{pipeline}.{step}.*", f"extraMem.{step}.*"]:
        for fp in glob.glob(os.path.join(stats_dir, pat)):
            try:
                os.remove(fp)
                print(f"Removed stat file: {fp}")
            except Exception as e:
                print(f"Could not remove {fp}: {e}")

# ── ASCII scatter plot (stdlib only) ─────────────────────────────────────────

def ascii_scatter(x_vals, y_vals, title="", xlabel="", ylabel="", width=88, height=22, deleted=None):
    """Draw a scatter plot in the terminal using only the standard library.
    Points are labeled 0-9 then a-z; collisions shown as '+'.
    deleted: set of indices to hide. Scale is always fixed to the full list so
             it does not change when dots are removed.
    """
    if deleted is None:
        deleted = set()

    Y_LBL  = 10   # chars reserved for y-axis tick labels
    plot_w = width - Y_LBL - 1   # 1 for the '|' separator
    plot_h = height

    # Scale is computed from ALL points (including deleted) so it never shifts.
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    if x_max == x_min: x_min -= 0.5; x_max += 0.5
    if y_max == y_min: y_min -= 0.5; y_max += 0.5

    xp = (x_max - x_min) * 0.05
    yp = (y_max - y_min) * 0.05
    x_min -= xp; x_max += xp
    y_min -= yp; y_max += yp

    grid  = [[' '] * plot_w for _ in range(plot_h)]
    chars = "0123456789abcdefghijklmnopqrstuvwxyz"

    for i, (x, y) in enumerate(zip(x_vals, y_vals)):
        if i in deleted:
            continue                     # skip deleted; index preserved for others
        gx = round((x - x_min) / (x_max - x_min) * (plot_w - 1))
        gy = round((y - y_min) / (y_max - y_min) * (plot_h - 1))
        gx = max(0, min(plot_w - 1, gx))
        gy = max(0, min(plot_h - 1, gy))
        gy = (plot_h - 1) - gy          # flip: y grows upward
        marker = chars[i] if i < len(chars) else '*'
        grid[gy][gx] = marker if grid[gy][gx] == ' ' else '+'

    if title:
        print(title.center(width))

    tick_rows = {0, plot_h // 2, plot_h - 1}
    for row_i, row in enumerate(grid):
        if row_i in tick_rows:
            y_frac = 1.0 - row_i / (plot_h - 1)
            y_val  = y_min + y_frac * (y_max - y_min)
            prefix = f"{y_val:{Y_LBL}.3f}|"
        else:
            prefix = ' ' * Y_LBL + '|'
        print(prefix + ''.join(row))

    # x axis
    print(' ' * Y_LBL + '+' + '-' * plot_w)

    x_lo  = f"{x_min:.3f}"
    x_mid = f"{(x_min + x_max) / 2:.3f}"
    x_hi  = f"{x_max:.3f}"
    gap1  = max(1, plot_w // 2 - len(x_lo) - len(x_mid) // 2)
    gap2  = max(1, plot_w - len(x_lo) - gap1 - len(x_mid) - len(x_hi))
    print(' ' * (Y_LBL + 1) + x_lo + ' ' * gap1 + x_mid + ' ' * gap2 + x_hi)

    if xlabel:
        print((' ' * (Y_LBL + 1) + xlabel).center(width))
    if ylabel:
        print(f"  y-axis: {ylabel}")

# ── Display ────────────────────────────────────────────────────────────────────

def show_plot_and_table(all_rows, label, deleted=None):
    """all_rows is the original full snapshot for this group (never trimmed).
    deleted is the set of indices that have been removed — they are hidden from
    the plot and table but their indices are preserved for all remaining points.
    """
    if deleted is None:
        deleted = set()

    if not all_rows:
        print(f"\n  (no records for '{label}')\n")
        return

    x_vals, y_vals, bad = [], [], []
    for i, row in enumerate(all_rows):
        try:
            x_vals.append(float(row["inputSize"]) / (1024 * 1024))  # bytes → GiB
            y_vals.append(float(row["memUsed"]) / 1024)             # MiB  → GiB
        except (ValueError, KeyError, TypeError):
            bad.append(i)
            x_vals.append(0.0)
            y_vals.append(0.0)

    remaining = len(all_rows) - len(deleted)
    if remaining == 0:
        print(f"\n  (all records for '{label}' have been deleted)\n")
        return

    ascii_scatter(
        x_vals, y_vals,
        title=f"Input Size (G) vs Memory (G) — {label}",
        xlabel="Input Size (G)",
        ylabel="Memory Usage (G)",
        deleted=deleted,
    )

    print(f"\n  Group: {label}   ({remaining} remaining of {len(all_rows)} total)")
    print(f"  {'Idx':>4}  {'Job ID':>12}  {'Input (G)':>12}  {'Mem (G)':>10}  {'Status':>10}")
    print("  " + "-" * 58)
    for i, (row, x, y) in enumerate(zip(all_rows, x_vals, y_vals)):
        if i in deleted:
            continue
        flag = " !" if i in bad else ""
        print(f"  {i:>4}  {str(row.get('jobID','')):>12}  {x:>12.4f}  {y:>10.4f}  {str(row.get('status','')):>10}{flag}")

# ── Index parsing ──────────────────────────────────────────────────────────────

def parse_indices(s, max_idx, exclude=None):
    """Parse '0', '1,3', '2-5', or combinations into a set of valid ints.
    exclude: set of indices already deleted (reported as invalid).
    """
    if exclude is None:
        exclude = set()
    indices = set()
    for part in s.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part:
            try:
                a, b = part.split("-", 1)
                indices.update(range(int(a), int(b) + 1))
            except ValueError:
                print(f"  Invalid range: {part}")
        else:
            try:
                indices.add(int(part))
            except ValueError:
                print(f"  Invalid index: {part}")
    valid = set()
    for i in indices:
        if i in exclude:
            print(f"  Index {i} already deleted — skipping.")
        elif 0 <= i < max_idx:
            valid.add(i)
        else:
            print(f"  Index {i} out of range — skipping.")
    return valid

# ── Main loop ──────────────────────────────────────────────────────────────────

def main():
    headers, rows = read_job_records()
    current_key   = None
    group_rows    = []     # full snapshot of rows for the selected group (never trimmed)
    deleted_indices = set()  # indices into group_rows that have been deleted

    while True:

        # ── Group selection ──────────────────────────────────────────────────
        if current_key is None:
            groups = unique_groups(rows)
            print("\n══════════════ SLURM Job Records — pipeline / step / ref ══════════════")
            for i, key in enumerate(groups):
                n = sum(1 for r in rows if group_key(r) == key)
                print(f"  [{i:>2}] {group_label(key):<50} {n:>4} rec")
            print("\n  Commands: <number> select group | s save | q quit")
            choice = input("  > ").strip().lower()

            if choice == "q":
                print("Bye.")
                break
            elif choice == "s":
                save_to_file(headers, rows)
                continue
            else:
                try:
                    current_key   = groups[int(choice)]
                    group_rows    = [r for r in rows if group_key(r) == current_key]
                    deleted_indices = set()
                except (ValueError, IndexError):
                    print("  Invalid choice — enter a number from the list.")
                    continue

        label = group_label(current_key)

        # ── Show plot and prompt for deletion ────────────────────────────────
        show_plot_and_table(group_rows, label, deleted=deleted_indices)

        print(
            "\n  Commands: <indices> delete (e.g. 0  or  1,3  or  2-5)"
            " | a delete all | b back | s save | q quit"
        )
        action = input("  > ").strip().lower()

        pipeline, step, _ref = current_key

        if action == "q":
            print("Bye.")
            break
        elif action == "b":
            current_key     = None
            group_rows      = []
            deleted_indices = set()
        elif action == "a":
            confirm = input(
                f"  Delete ALL {len(group_rows) - len(deleted_indices)} remaining"
                f" records for '{label}'? [y/N] "
            ).strip().lower()
            if confirm == "y":
                all_indices = set(range(len(group_rows))) - deleted_indices
                job_ids = {group_rows[i]["jobID"] for i in all_indices}
                deleted_indices |= all_indices
                rows = [r for r in rows if r["jobID"] not in job_ids]
                write_tmp(headers, rows)
                delete_stats_files(pipeline, step)
                print(f"  Deleted all {len(job_ids)} records for '{label}'.")
            else:
                print("  Cancelled.")
        elif action == "s":
            save_to_file(headers, rows)
        elif action:
            indices = parse_indices(action, len(group_rows), exclude=deleted_indices)
            if indices:
                job_ids = {group_rows[i]["jobID"] for i in indices}
                deleted_indices |= indices
                rows = [r for r in rows if r["jobID"] not in job_ids]
                write_tmp(headers, rows)
                delete_stats_files(pipeline, step)
                print(
                    f"  Deleted {len(job_ids)} job(s): "
                    + ", ".join(sorted(job_ids))
                )
            else:
                print("  No valid indices — try again.")
        # empty input → just redraw

if __name__ == "__main__":
    main()
