# Demo Nextflow Pipeline with nf-diff CI

A small but realistic FASTQ QC/trimming Nextflow pipeline configured with [nf-diff](https://github.com/mribeirodantas/nf-diff) to automatically test and diff runs on Pull Requests against the baseline run on `main`.

## Pipeline Structure

- `data/`: Three synthetic single-end FASTQ samples (`sampleA`, `sampleB`, `sampleC`) with reads of varying length and GC content.
- `main.nf`: A 5-process workflow:
  - `RAW_QC` / `TRIMMED_QC`: compute read count, total bases, mean read length, and GC% from a FASTQ file.
  - `TRIM_FILTER`: drops reads shorter than `params.min_length` and trims `params.trim_bases` off the 3' end of the rest.
  - `SAMPLE_REPORT`: combines the raw and trimmed QC stats per sample (read retention %, length/GC before and after).
  - `AGGREGATE_REPORT`: collects all per-sample rows into one `summary_report.tsv`.
- `nextflow.config`: Defines the sample sheet (`params.samples`), `min_length`, `trim_bases`, `outdir`, and a `test` profile.
- `.github/workflows/update-baseline.yml`: Runs on pushes to `main` and caches the `.nextflow` run cache and baseline output results.
- `.github/workflows/pr-diff.yml`: Runs on `pull_request`, restores the baseline `.nextflow` cache from `main`, runs the PR pipeline, and posts a rich Markdown diff report as a PR comment using `nf-diff`.

## Testing nf-diff in Action

1. Push this repository to GitHub as your `main` branch.
2. The `Update Baseline Run` workflow will trigger and populate the initial baseline cache.
3. Create a feature branch:
   ```bash
   git checkout -b feature/adjust-trim-params
   ```
4. Introduce a realistic change, such as:
   - Raising `min_length` or `trim_bases` in `nextflow.config` (changes how many reads survive filtering and the resulting mean length/GC% per sample).
   - Adding a new sample to `params.samples`.
   - Tweaking the QC/trim `awk` logic in `main.nf`.
5. Push the branch and open a Pull Request.
6. Observe the `PR nf-diff Check` workflow running Nextflow once and posting the interactive comparison comment on your PR, with per-sample differences in read counts, retention %, mean length, and GC%!
