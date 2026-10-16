#!/usr/bin/env python3
"""Run every code check of this repository and compare it with the stored output.

For every folder ``code/NN_slug/`` this script

1. collects the Python scripts of the folder,
2. runs each one from inside its folder (``cd NN_slug && python3 script.py``),
3. compares what it prints on stdout with ``NN_slug/expected_output/<script>.txt``,
4. prints a pass/fail table at the end.

Special cases
    * ``09_schwinger_keldysh/kexp.py`` takes the expansion order N as its argument.  It is run once for
      every stored file ``expected_output/kexp_N<N>.txt`` (N = 2, 4, 6, 10), as in that folder's
      ``run_all.sh``.
    * Helper modules without a stored output (``engine.py`` and ``rows.py`` in ``09_schwinger_keldysh``)
      are executed as well; they must exit cleanly and print nothing.
    * A stored output whose script is missing is reported as ``ERROR`` instead of being skipped.
    * Several scripts import a sibling module (for example ``check_global.py``), so the scripts are run
      with their own folder as working directory and *without* Python's isolated mode (``-I``).

Comparison
    Line endings are normalized, trailing whitespace on each line and trailing blank lines are ignored,
    and pure timing lines (for example ``elapsed: 12.3 s``) are dropped on both sides.  Nothing else is
    normalized: every number, sign and symbolic expression must agree exactly.  The stored outputs were
    produced with Python 3.13.1, SymPy 1.14.0 and mpmath 1.3.0; another SymPy version may print an
    equivalent expression in a different form, which shows up here as a mismatch (see the diff).

Usage (from the repository root or from ``code/``)
    python3 code/run_all.py               # everything: 38 runs, about 10 minutes on a laptop
    python3 code/run_all.py --quick       # skip the two scripts that take longer than 60 s (about 6 minutes)
    python3 code/run_all.py --folder 09   # one folder (can be repeated: --folder 02 --folder 08)
    python3 code/run_all.py --list        # show what would be run, with estimated runtimes
    python3 code/run_all.py -j 4          # run up to 4 scripts in parallel

Options ``--keep DIR`` (save every actual output to DIR) and ``--timeout SECONDS`` are also available.
The exit status is 0 if every executed script matches its stored output, and 1 otherwise.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import difflib
import os
import re
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

CODE_DIR = Path(__file__).resolve().parent

# Approximate runtimes in seconds, measured in a full run on an Apple-silicon Mac (Python 3.13.1, SymPy 1.14.0).
# Used for --list and for --quick.  A script that is not listed (e.g. a newly added one) is assumed to be fast.
RUNTIME_ESTIMATE: dict[str, float] = {
    "01_feynman_rules/check_eq6.py": 22,
    "01_feynman_rules/check_vertices.py": 10,
    "02_four_point/check_41_44.py": 1,
    "02_four_point/check_4pt.py": 23,
    "02_four_point/check_coincidence.py": 2,
    "03_three_point_setup/check_45_46.py": 171,
    "03_three_point_setup/check_47_49.py": 2,
    "03_three_point_setup/check_52.py": 1,
    "04_three_point_local_log/check_54_59.py": 1,
    "04_three_point_local_log/check_60_67.py": 8,
    "05_normal_terms_tableIII/check_68_80.py": 17,
    "05_normal_terms_tableIII/check_81_93.py": 19,
    "06_tables_IV_V/check_94_125.py": 92,
    "06_tables_IV_V/check_global.py": 13,
    "06_tables_IV_V/check_table4.py": 2,
    "06_tables_IV_V/check_tables45.py": 20,
    "06_tables_IV_V/check_tables45b.py": 30,
    "06_tables_IV_V/diag_leading_pair.py": 2,
    "06_tables_IV_V/superseded_test_E.py": 15,
    "07_desitter_counterterms/check_126_135.py": 9,
    "07_desitter_counterterms/check_136_145.py": 5,
    "08_renormalization_tableX/check_146_IX.py": 2,
    "08_renormalization_tableX/check_160_167.py": 2,
    "08_renormalization_tableX/check_global.py": 13,
    "08_renormalization_tableX/check_tableX_global.py": 17,
    "09_schwinger_keldysh/analysis2.py": 10,
    "09_schwinger_keldysh/analysis3.py": 11,
    "09_schwinger_keldysh/engine.py": 0.2,
    "09_schwinger_keldysh/extract.py": 35,
    "09_schwinger_keldysh/kexp.py 2": 2,
    "09_schwinger_keldysh/kexp.py 4": 5,
    "09_schwinger_keldysh/kexp.py 6": 10,
    "09_schwinger_keldysh/kexp.py 10": 25,
    "09_schwinger_keldysh/numcheck.py": 4,
    "09_schwinger_keldysh/numcheck_k.py": 5,
    "09_schwinger_keldysh/rows.py": 0.2,
    "09_schwinger_keldysh/table5.py": 1,
    "09_schwinger_keldysh/tables234.py": 1,
}
QUICK_LIMIT = 60.0  # seconds; --quick skips every job whose estimate exceeds this

# Lines that contain nothing but a timing report are ignored in the comparison.
TIMING_LINE = re.compile(
    r"^\s*[\[(]?\s*(runtime|elapsed|time|took|wall(\s*time)?|cpu(\s*time)?)\b[^=]*?[:=]?\s*"
    r"[0-9]+(\.[0-9]+)?\s*(s|sec|secs|seconds|ms|min)\.?\s*[\])]?\s*$",
    re.IGNORECASE,
)

PARAM_FILE = re.compile(r"^(?P<stem>.+)_N(?P<n>\d+)\.txt$")  # e.g. kexp_N10.txt -> kexp.py 10


@dataclass
class Job:
    folder: str
    script: str
    args: list[str] = field(default_factory=list)
    expected: Path | None = None  # None: helper module, must print nothing

    @property
    def key(self) -> str:
        return f"{self.folder}/{self.script}" + ("".join(" " + a for a in self.args))

    @property
    def label(self) -> str:
        return self.script + ("".join(" " + a for a in self.args))

    @property
    def estimate(self) -> float:
        return RUNTIME_ESTIMATE.get(self.key, 5.0)


@dataclass
class Result:
    job: Job
    status: str  # PASS, FAIL, ERROR, TIMEOUT, SKIP
    seconds: float = 0.0
    detail: str = ""


def normalize(text: str) -> list[str]:
    lines = [ln.rstrip() for ln in text.replace("\r\n", "\n").replace("\r", "\n").split("\n")]
    lines = [ln for ln in lines if not TIMING_LINE.match(ln)]
    while lines and not lines[-1]:
        lines.pop()
    return lines


def collect_jobs(folder: Path) -> list[Job]:
    exp_dir = folder / "expected_output"
    expected = {p.name for p in exp_dir.glob("*.txt")} if exp_dir.is_dir() else set()
    jobs: list[Job] = []
    used: set[str] = set()
    for script in sorted(folder.glob("*.py")):
        stem = script.stem
        if f"{stem}.txt" in expected:
            jobs.append(Job(folder.name, script.name, [], exp_dir / f"{stem}.txt"))
            used.add(f"{stem}.txt")
            continue
        params = sorted(
            (int(m.group("n")), name)
            for name in expected
            if (m := PARAM_FILE.match(name)) and m.group("stem") == stem
        )
        if params:
            for n, name in params:
                jobs.append(Job(folder.name, script.name, [str(n)], exp_dir / name))
                used.add(name)
        else:
            jobs.append(Job(folder.name, script.name, [], None))  # helper module
    # A stored output without its script is not silently ignored: it becomes a job for the missing
    # script, which then fails with ERROR ("can't open file ...").
    for name in sorted(expected - used):
        m = PARAM_FILE.match(name)
        stem, args = (m.group("stem"), [m.group("n")]) if m else (name[: -len(".txt")], [])
        jobs.append(Job(folder.name, f"{stem}.py", args, exp_dir / name))
    return jobs


def run_job(job: Job, timeout: float, keep: Path | None) -> Result:
    cwd = CODE_DIR / job.folder
    env = dict(os.environ)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    env["PYTHONIOENCODING"] = "utf-8"
    t0 = time.perf_counter()
    try:
        proc = subprocess.run(
            [sys.executable, job.script, *job.args],
            cwd=cwd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return Result(job, "TIMEOUT", time.perf_counter() - t0, f"no result after {timeout:.0f} s")
    dt = time.perf_counter() - t0
    out = proc.stdout.decode("utf-8", errors="replace")
    err = proc.stderr.decode("utf-8", errors="replace")

    if keep is not None:
        target = keep / job.folder
        target.mkdir(parents=True, exist_ok=True)
        name = job.expected.name if job.expected else f"{Path(job.script).stem}.txt"
        (target / name).write_text(out, encoding="utf-8")

    if proc.returncode != 0:
        tail = "\n".join(err.strip().splitlines()[-15:])
        return Result(job, "ERROR", dt, f"exit status {proc.returncode}\n{tail}")

    want_text = job.expected.read_text(encoding="utf-8") if job.expected else ""
    got, want = normalize(out), normalize(want_text)
    if got == want:
        return Result(job, "PASS", dt)

    diff = list(
        difflib.unified_diff(
            want, got,
            fromfile=f"expected_output/{job.expected.name}" if job.expected else "(no output expected)",
            tofile="actual stdout",
            lineterm="", n=1,
        )
    )
    shown = "\n".join(diff[:40]) + ("\n  ..." if len(diff) > 40 else "")
    if err.strip():
        shown += "\n  stderr: " + "\n  ".join(err.strip().splitlines()[-5:])
    return Result(job, "FAIL", dt, shown)


def select_folders(requested: list[str] | None) -> list[Path]:
    folders = sorted(p for p in CODE_DIR.iterdir() if p.is_dir() and re.match(r"^\d\d_", p.name))
    if not requested:
        return folders
    chosen = []
    for req in requested:
        req = req.strip().rstrip("/")
        hits = [p for p in folders if p.name == req or p.name.startswith(req.zfill(2) + "_")]
        if not hits:
            sys.exit(f"run_all.py: no code folder matches '{req}'. "
                     f"Available: {', '.join(p.name for p in folders)}")
        chosen.extend(h for h in hits if h not in chosen)
    return chosen


def fmt_time(sec: float) -> str:
    return f"{sec:6.1f} s" if sec < 600 else f"{sec / 60:5.1f} min"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quick", action="store_true",
                    help=f"skip scripts estimated to take longer than {QUICK_LIMIT:.0f} s")
    ap.add_argument("--folder", action="append", metavar="NN",
                    help="run only this folder, e.g. 09 or 09_schwinger_keldysh (repeatable)")
    ap.add_argument("--list", action="store_true", help="list the jobs and exit")
    ap.add_argument("-j", "--jobs", type=int, default=1, help="number of scripts to run in parallel")
    ap.add_argument("--timeout", type=float, default=3600, help="per-script timeout in seconds")
    ap.add_argument("--keep", type=Path, metavar="DIR", help="save every actual output under DIR")
    opts = ap.parse_args()

    jobs = [j for f in select_folders(opts.folder) for j in collect_jobs(f)]

    if opts.list:
        for j in jobs:
            what = j.expected.name if j.expected else "(helper module, no output)"
            mark = "  [long, skipped by --quick]" if j.estimate > QUICK_LIMIT else ""
            print(f"{j.folder:28s} {j.label:24s} ~{max(j.estimate, 1):4.0f} s  -> {what}{mark}")
        return 0

    to_run = [j for j in jobs if not (opts.quick and j.estimate > QUICK_LIMIT)]
    skipped = [Result(j, "SKIP", 0.0, f"--quick (estimated ~{j.estimate:.0f} s)")
               for j in jobs if j not in to_run]
    est = sum(j.estimate for j in to_run)
    est_txt = f"~{est:.0f} s" if est < 90 else f"~{est / 60:.0f} min"
    print(f"Running {len(to_run)} script(s) from {CODE_DIR} with {sys.executable} "
          f"(Python {sys.version.split()[0]}); estimated total {est_txt}"
          + (f" sequentially, less with -j {opts.jobs}" if opts.jobs > 1 else "") + ".")
    if skipped:
        print(f"Skipping {len(skipped)} long script(s) because of --quick.")
    print(flush=True)

    keep = opts.keep.resolve() if opts.keep else None
    results: dict[str, Result] = {r.job.key: r for r in skipped}
    t_start = time.perf_counter()

    def report(r: Result) -> None:
        print(f"  {r.status:7s} {fmt_time(r.seconds)}  {r.job.folder}/{r.job.label}", flush=True)

    if opts.jobs > 1:
        with cf.ThreadPoolExecutor(max_workers=opts.jobs) as pool:
            futs = {pool.submit(run_job, j, opts.timeout, keep): j for j in to_run}
            for fut in cf.as_completed(futs):
                r = fut.result()
                results[r.job.key] = r
                report(r)
    else:
        for j in to_run:
            r = run_job(j, opts.timeout, keep)
            results[r.job.key] = r
            report(r)
    total = time.perf_counter() - t_start

    # Summary table in folder/script order.
    ordered = [results[j.key] for j in jobs]
    w_f = max(len("folder"), *(len(r.job.folder) for r in ordered))
    w_s = max(len("script"), *(len(r.job.label) for r in ordered))
    w_e = max(len("compared with"), *(len(r.job.expected.name if r.job.expected else "(no output)")
                                      for r in ordered))
    line = f"+-{'-' * w_f}-+-{'-' * w_s}-+-{'-' * w_e}-+---------+----------+"
    print("\n" + line)
    print(f"| {'folder':{w_f}} | {'script':{w_s}} | {'compared with':{w_e}} | result  |     time |")
    print(line)
    for r in ordered:
        exp = r.job.expected.name if r.job.expected else "(no output)"
        t = fmt_time(r.seconds) if r.status != "SKIP" else "       -"
        print(f"| {r.job.folder:{w_f}} | {r.job.label:{w_s}} | {exp:{w_e}} | {r.status:7s} | {t:>8s} |")
    print(line)

    counts = {s: sum(r.status == s for r in ordered) for s in ("PASS", "FAIL", "ERROR", "TIMEOUT", "SKIP")}
    print("  ".join(f"{k}: {v}" for k, v in counts.items() if v or k == "PASS")
          + f"   (wall time {fmt_time(total).strip()})")

    bad = [r for r in ordered if r.status in ("FAIL", "ERROR", "TIMEOUT")]
    for r in bad:
        print(f"\n--- {r.status}: {r.job.folder}/{r.job.label}\n{r.detail}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
