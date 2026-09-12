# Demo Nextflow Pipeline with nf-diff CI

A simple Nextflow pipeline repository configured with [nf-diff](https://github.com/mribeirodantas/nf-diff) to automatically test and diff runs on Pull Requests against the baseline run on `main`.

## Pipeline Structure

- `main.nf`: Processes text input channels across `SPLIT_WORDS`, `COUNT_WORDS`, and `UPPERCASE` processes.
- `nextflow.config`: Contains default parameters and a `test` profile.
- `.github/workflows/update-baseline.yml`: Runs on pushes to `main` and caches the `.nextflow` run cache and baseline output results.
- `.github/workflows/pr-diff.yml`: Runs on `pull_request`, restores the baseline `.nextflow` cache from `main`, runs the PR pipeline, and posts a rich Markdown diff report as a PR comment using `nf-diff`.

## Testing nf-diff in Action

1. Push this repository to GitHub as your `main` branch.
2. The `Update Baseline Run` workflow will trigger and populate the initial baseline cache.
3. Create a feature branch:
   ```bash
   git checkout -b feature/test-diff
   ```
4. Introduce a change, such as:
   - Modifying a process script in `main.nf`.
   - Changing default parameters in `nextflow.config`.
   - Adding a new process or modifying channel transformations.
5. Push the branch and open a Pull Request.
6. Observe the `PR nf-diff Check` workflow running Nextflow once and posting the interactive comparison comment on your PR!
