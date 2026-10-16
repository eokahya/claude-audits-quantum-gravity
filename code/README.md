# Code checks

This folder holds the independent computer-algebra and numerical checks of the re-derivation. There is one folder per stage
of the calculation, numbered like the [docs pages](../docs/) and the chapters of the [notes](../notes/notes.pdf).
The scripts **check** derivations that were first done analytically. They do not replace those derivations.
Finding IDs such as S2 or C5 below refer to the [complete list of findings](../docs/10-complete-list-of-findings.md).

Every script runs on its own from its folder, prints its results, and has a stored output in `expected_output/`. That stored
output is identical to the output recorded when the calculation was originally carried out.

## Requirements

- Python 3.10 or later
- `sympy` and `mpmath` (`pip install -r requirements.txt`); no script needs `numpy` or `scipy`
- no network access, no other files; the scripts read no data and write nothing to disk

The stored outputs were produced with **Python 3.13.1, SymPy 1.14.0 and mpmath 1.3.0** on macOS. With other SymPy versions the
mathematics is the same, but an expression may be printed in a different (equivalent) form. For an exact comparison use
`pip install sympy==1.14.0 mpmath==1.3.0`.

## Run everything

From the repository root (inside `code/`, drop the `code/` prefix):

```bash
python3 code/run_all.py              # all 38 runs, about 10 minutes on an Apple-silicon laptop
python3 code/run_all.py --quick      # skip the two scripts that take more than 60 s (about 6 minutes)
python3 code/run_all.py --folder 06  # one folder; repeat the option for several (--folder 02 --folder 08)
python3 code/run_all.py --list       # show the runs and their estimated runtimes without running anything
python3 code/run_all.py -j 4         # run up to four scripts in parallel
```

| Option | Meaning |
|---|---|
| `--quick` | Skip scripts whose estimated runtime is above 60 s: `03_three_point_setup/check_45_46.py` (~3 min) and `06_tables_IV_V/check_94_125.py` (~1.5 min). |
| `--folder NN` | Run only folder `NN` (`06` or `06_tables_IV_V`). Can be repeated. |
| `--list` | Print every run, the file it is compared with, and its estimated runtime. |
| `-j N` | Run N scripts at the same time. The summary table keeps the folder order. |
| `--keep DIR` | Also save every actual output under `DIR/<folder>/`. Useful for inspecting a mismatch. |
| `--timeout S` | Per-script time limit in seconds (default 3600). |

The script prints one line per run as it finishes, then a pass/fail table, and then a diff for every run that failed. Its exit
status is 0 only if every run passed.

**Last full run** (16 October 2026, Apple-silicon Mac, Python 3.13.1, SymPy 1.14.0, mpmath 1.3.0): `PASS: 38` with no
failures, in 10.1 minutes of wall time.

## How the outputs are compared

`run_all.py` does the following for every folder `NN_slug/`:

1. Every `*.py` file whose stem has a stored output `expected_output/<stem>.txt` is run as
   `cd NN_slug && python3 <stem>.py`, with the same interpreter that runs `run_all.py`.
2. A script with stored outputs of the form `expected_output/<stem>_N<n>.txt` takes an argument. It is run once for each `n`.
   This applies only to `09_schwinger_keldysh/kexp.py`, which is run as `python3 kexp.py 2`, `4`, `6` and `10`, as in that
   folder's `run_all.sh`.
3. A script without any stored output is a helper module that other scripts import (`engine.py` and `rows.py` in
   `09_schwinger_keldysh`). It is run as well, and it must exit cleanly and print nothing.
   Conversely, a stored output whose script is missing is reported as `ERROR`, so no stored output is skipped silently.
4. The scripts are **not** run in Python's isolated mode (`-I`), because several of them import a sibling module from their
   own folder: `check_global.py` in `06` and `08`, and `engine.py`, `rows.py`, `kexp.py`, `tables234.py` and `analysis2.py` in `09`.

Before comparing, only these differences are normalized: line endings, trailing whitespace on each line, trailing blank lines,
and lines that contain nothing but a timing report (the current scripts print none). Every number, sign and printed
expression must agree exactly. A script that exits with an error counts as `ERROR`. Text on stderr is not compared, but
it is shown when a run fails.

Note: `09_schwinger_keldysh/run_all.sh` *regenerates* the stored outputs of that folder, so it overwrites them.
To *compare* with the stored outputs, use `run_all.py`.

## The folders

| Folder | Stage (docs page) | What it checks | Scripts | Runtime |
|---|---|---|---|---|
| [`01_feynman_rules`](01_feynman_rules/) | [1. Feynman rules](../docs/01-feynman-rules.md) | (I.4)–(I.6) against the exact inverse metric and determinant. Tables I–II from the exact metric alone, including the sign of $`U_2`$ (C5). | `check_eq6.py`, `check_vertices.py` | ~30 s |
| [`02_four_point`](02_four_point/) | [2. 4-point contribution](../docs/02-four-point.md) | Lemma A and the overall sign of (I.33) (S1). Coincidence limits (I.34)–(I.40), including the sign of the finite part of (I.40) (B1). The assembly (I.41)–(I.44). | `check_4pt.py`, `check_coincidence.py`, `check_41_44.py` | ~25 s |
| [`03_three_point_setup`](03_three_point_setup/) | [3. 3-point setup](../docs/03-three-point-setup.md) | (I.45) → (I.46) and the derivative rule. The nine terms (I.47)–(I.49) equal $`-1`$ times the printed ones (S1). The local term of (I.52) and the first coefficient of (I.14) (C4). | `check_45_46.py` (long), `check_47_49.py`, `check_52.py` | ~3 min |
| [`04_three_point_local_log`](04_three_point_local_log/) | [4. local and log terms](../docs/04-three-point-local-and-log.md) | Local terms (I.54)–(I.59), including $`\tfrac72\to\tfrac52`$ (B2). Logarithm terms (I.60)–(I.67) and their cancellation against (I.44). | `check_54_59.py`, `check_60_67.py` | ~10 s |
| [`05_normal_terms_tableIII`](05_normal_terms_tableIII/) | [5. normal terms, Table III](../docs/05-normal-terms-table-III.md) | The derivative-extraction and time-derivative identities (I.68)–(I.93). Table III, and the missing $`H^2`$ in (I.89) (C8). | `check_68_80.py`, `check_81_93.py` | ~40 s |
| [`06_tables_IV_V`](06_tables_IV_V/) | [6. Tables IV–V](../docs/06-tables-IV-V.md) | (I.94)–(I.125). The appendix totals against a brute-force direct sum. Tables IV + V: the sign of $`\gamma_2`$ (C15) and the missing $`72\gamma_2`$ (S3). | `check_94_125.py` (long), `check_global.py`, `check_tables45.py`, `check_tables45b.py`, `check_table4.py`, `diag_leading_pair.py`, `superseded_test_E.py` | ~3 min |
| [`07_desitter_counterterms`](07_desitter_counterterms/) | [7. symmetry and counterterms](../docs/07-de-sitter-symmetry-and-counterterms.md) | Killing vectors of de Sitter in conformal coordinates. Which of them the gauge breaks. The basis of counterterms (I.136)–(I.145). | `check_126_135.py`, `check_136_145.py` | ~15 s |
| [`08_renormalization_tableX`](08_renormalization_tableX/) | [8. renormalized self-mass](../docs/08-renormalized-self-mass-table-X.md) | (I.146)–(I.159), Tables VI–IX, $`\beta_1,\beta_2,\beta_3`$ (B3, B4). The assembly of Table X, and the global 30-digit pointwise test, which passes only with both corrections (S2 and S3). | `check_146_IX.py`, `check_160_167.py`, `check_tableX_global.py`, `check_global.py` (helper with a self-test) | ~40 s |
| [`09_schwinger_keldysh`](09_schwinger_keldysh/) | [9. Schwinger–Keldysh, late times](../docs/09-schwinger-keldysh-late-time.md) | Paper II Tables 2–5. The full right-hand side beyond (II.44) (S4), expanded exactly through $`k^{10}`$. The late-time coefficients. | `tables234.py`, `numcheck.py`, `table5.py`, `kexp.py N`, `numcheck_k.py`, `analysis2.py`, `analysis3.py`, `extract.py`; helpers `engine.py`, `rows.py` | ~2 min |

Runtimes were measured on an Apple-silicon Mac. Each folder's README lists every script, the equations it checks, and how to read
its output line by line.

## Conventions used in the outputs

- **Sign.** Unless a script says otherwise, it works in the *paper's convention*, the sign of the printed self-mass-squared.
  The *definition sign* fixed by Paper I's own Eq. (I.168) is $`-1`$ times this. Scripts that print both say so.
- **Random points.** The pointwise numerical checks use fixed seeds or fixed rational points, so their output is reproducible.
  They are checks at sample points (two or three spacelike point pairs, 30 digits), not symbolic identities.
- **$`\kappa^2`$ in folder 09.** The docs and notes write $`\Phi=u+\kappa^2\Phi_1`$, where the source $`\mathcal R_1`$ of
  $`\Phi_1`$ is built from the one-loop self-mass with $`\kappa^2`$ stripped off, and quote Table 5 entries and
  $`\mathcal R_1`$ in units of $`H^4u(0,k)/16\pi^2`$. The scripts use the same units (see the docstring of `table5.py`).
  The original working report wrote an extra $`\kappa^2`$ in this unit; the single $`\kappa^2`$ of
  $`\Phi=u+\kappa^2\Phi_1`$ must not be counted a second time. (The header $`\kappa^2H^4/(4\pi)^4`$ in `rows.py` is the
  unit of the Table X / Table 1 coefficients, not of the source.)
- **"Expected" means "recorded".** A stored output is what the script printed in the original calculation; it is not a
  target. Some outputs deliberately show a `False`: for example, when the printed version of an equation is tested and fails.
  Each folder's README explains every such line.

## Provenance

The scripts are the ones that were run during the 18 tasks of the original calculation. For Tasks 6–18 they were taken from
the working directory. For Tasks 1–5 they were reconstructed from the command logs. Each folder README lists the changes
made while cleaning them up: a docstring at the top of each script, removal of absolute paths and unused code, cosmetic
reformatting, and a deterministic sort key in one script. In folder `09`, the pickle files that used to pass results between
scripts are gone (each script now recomputes what it needs), and `kexp.py` is wrapped in a class. None of these changes
touches the mathematics. Every stored output was regenerated from the cleaned script and compared with the original log, and
they agree.
