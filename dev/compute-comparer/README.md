# compute-comparer

Proves that the `.pgx` produced by the **Shiny upload flow** (omicsplayground `master`
+ playbase `main`) matches the one produced by the **script path** (omicsplayground
`edgy` + playbase `edgy`, using `pgx.createPGX(preprocess=)`), across a matrix of
preprocessing settings.

**Results so far:** the seam is exact. `pgx.preprocess()` reproduces the app's `countsX`
bit-for-bit across every normalization method, every NA-handling path, and outlier
removal — including the row- and column-reducing cases (7439→5531 features,
18→17 samples), verified through real browser sessions. See `runs/report.csv`.

Two pre-existing **app bugs** turned up along the way, written up in
[`../../docs/FINDINGS.md`](../../docs/FINDINGS.md): the batch-correction method is
silently reset to ComBat, and every wizard lock/unlock branch is dead code.

Background on why the two ever diverged, and the full mechanics:

- [`docs/app-automation.md`](docs/app-automation.md) — driving the upload wizard with Playwright
- [`docs/script-pipeline.md`](docs/script-pipeline.md) — the `params.RData` contract and the known implementation deltas

## The idea in one paragraph

The app does not compute the pgx itself. It writes `params.RData` into a temp dir and
shells out to `bin/pgxcreate_op.R`. That file is therefore the exact app↔script seam:
capture it, and you can reproduce the app's pgx offline. So Playwright only drives the
browser far enough to make `params.RData` appear (~2 min), and both pipelines are then
executed by us under fixed seeds and pinned library paths — instead of one of them
coming from an uncontrolled 30-60 min background job.

## Layout

```
drive_app.mjs      Path A: browser -> params.RData        (playbase main, master tree)
run_script.R       Path B: fixture CSVs + scenario -> params.RData   (playbase edgy)
run_pgx.R          executes either params.RData: createPGX only, or the full script
compare.R          the three gates; appends to runs/report.csv
scenarios.json     the settings matrix
libs/pb-main       playbase@main 4.1.3  (no pgx.preprocess)     -- gitignored
libs/pb-edgy       playbase@edgy 4.1.02 (has pgx.preprocess)    -- gitignored
runs/              per-scenario artifacts + report.csv          -- gitignored
```

## Gates

| Gate | Compares | Cost | Deterministic |
|---|---|---|---|
| 1 | `params.RData` shared fields, **and** `app countsX == pgx.preprocess(raw, opts)$X` | seconds | yes |
| 2 | `pgx.createPGX()` output: counts, X, samples, contrasts, genes | minutes | yes |
| 3 | the full `.pgx` | 30-60 min | no |

Gate 1 is the sharp one. Master's `params` has no `preprocess` field, so gate 1 cannot
diff the mapping directly — instead it re-runs `pgx.preprocess()` on the raw counts with
the scenario's settings and checks the result against the `countsX` the **real browser
session** produced. That isolates the seam completely: no `createPGX`, no compute, no RNG.

## Running

```bash
make check                  # confirm both playbase trees resolve as expected
make app.start              # launch the app under playbase main (port 3838)

make gate1 SC=baseline      # = a + b + g1
make gate2 SC=baseline      # = create.a + create.b + g2

make app.stop
```

Individual steps: `make a` / `make b` (capture), `make create.a` / `make create.b`,
`make full.a` / `make full.b`, `make g1` / `g2` / `g3`. `SC` selects the scenario,
`SEED` the RNG seed (default 42).

One-time setup (already done in this worktree): `make setup`.

## Prerequisites

- git worktrees at `../../../omicsplayground-edgy` (branch `edgy`) and
  `../../../playbase-main` (branch `main`)
- `etc/OPTIONS` with `AUTHENTICATION = none` (the driver has no login leg)
- `make sass` run once in the omicsplayground root so `components/00SourceAll.R` exists

## Reading the results

`runs/report.csv` gets one row per check: `scenario, gate, check, status, detail`.

- `PASS` / `FAIL` — the verdict.
- `INFO` — a gate-3 slot known to be stochastic (t-SNE, UMAP, clustering, WGCNA).
  `pgx.computePGX` seeds nothing and forked workers ignore the top-level seed, so these
  differ run to run. **Establish the noise floor first** by running the same scenario
  twice on Path A and diffing A against A — otherwise there is no way to tell divergence
  from jitter.
- `SKIP` / `WARN` — usually a missing input or an expected field-inventory difference.

## Caveat on gate 3

playbase `main` (4.1.3) and `edgy` (4.1.02) differ by far more than the preprocess work,
so a gate-3 difference is not automatically a preprocess bug. Gates 1 and 2 are what
isolate the seam; gate-3 differences need attribution before they mean anything.
