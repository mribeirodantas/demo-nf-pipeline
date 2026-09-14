# Demo Nextflow Pipeline with nf-diff CI

A small but realistic FASTQ QC/trimming Nextflow pipeline configured with [nf-diff](https://github.com/mribeirodantas/nf-diff) to automatically test and diff runs on Pull Requests against the baseline run on `main`.

## Pipeline Structure

- `data/`: Three synthetic single-end FASTQ samples (`sampleA`, `sampleB`, `sampleC`) with reads of varying length and GC content.
- `main.nf`: A 5-process workflow:
  - `RAW_QC`: computes read count, total bases, mean read length, and GC% from the raw FASTQ.
  - `LENGTH_HISTOGRAM`: buckets raw reads into short (<30bp) / medium (30-59bp) / long (>=60bp) counts.
  - `TRIM_AND_QC`: trims `params.trim_bases5` off the 5' end and `params.trim_bases` off the 3' end, drops reads shorter than `params.min_length`, and computes QC stats on the trimmed reads in one step.
  - `SAMPLE_REPORT`: combines raw QC, trimmed QC, and the length histogram per sample (read retention %, length/GC before and after, read-length bucket counts).
  - `AGGREGATE_REPORT`: collects all per-sample rows into one `summary_report.tsv`.
- `nextflow.config`: Defines the sample sheet (`params.samples`), `min_length`, `trim_bases`, `trim_bases5`, `outdir`, and a `test` profile.
- `.github/workflows/update-baseline.yml`: Runs on pushes to `main` and caches the `.nextflow` run cache and baseline output results.
- `.github/workflows/pr-diff.yml`: Runs on `pull_request`, restores the baseline `.nextflow` cache from `main`, runs the PR pipeline, and posts a rich Markdown diff report as a PR comment using `nf-diff`.

## Testing nf-diff in Action

1. Push this repository to GitHub as your `main` branch.
2. The `Update Baseline Run` workflow will trigger and populate the initial baseline cache.
3. Create a feature branch:
   ```bash
   git checkout -b feature/adjust-trim-params
   ```
4. Introduce a structural change to the pipeline itself, such as:
   - Adding a new process (e.g. a length-histogram or adapter-scan step) and wiring its output into the reports.
   - Removing or merging processes (e.g. combining separate trim and QC steps into one).
   - Rewriting a process's `awk`/script logic (e.g. changing trimming from 3'-only to 5'+3', or changing how GC%/length are computed).
5. Push the branch and open a Pull Request.
6. Observe the `PR nf-diff Check` workflow running Nextflow once and posting the interactive comparison comment on your PR, showing added/removed processes, new or missing output files, and per-sample differences in read counts, retention %, mean length, and GC%!
