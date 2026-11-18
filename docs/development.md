# Development and publishing

For running analyses from a published build, start with the [README](../README.md).

## Build from source

Use Viash 0.9.7, Nextflow 25.10.6, Java 17 through 24 and Docker.
Clone the source into a separate checkout from any published build:

```bash
git clone https://github.com/sunneyxie-lab/cpl-sc-foodie.git /absolute/path/to/cpl-sc-foodie-source
repo=/absolute/path/to/cpl-sc-foodie-source
(cd "$repo" && viash ns build --parallel --setup cachedbuild)
```

Run the generated workflows under `target/nextflow/workflows/`, rather than the
source files under `src/`. Keep the complete generated tree and its dependencies,
including the configuration files beside each `main.nf`.

Source branches, even ones named like versions, do not contain generated
`target/`. Set `genome` to your prepared Bismark genome and follow the analysis
commands in the [README](../README.md#process-multiple-samples).

## CI builds and releases

### Workflows and outputs

CI uses Viash 0.9.7, Nextflow 25.10.6, Java 17 and Docker on standard
`ubuntu-latest` runners.

| Workflow | Trigger and responsibility |
| --- | --- |
| `viash-test.yml` | Pull requests, pushes to `main`, or reusable calls: all local production components' synthetic tests. PRs are read-only and never publish. |
| `main-build.yml` | Pushes to `main`, or reusable integration/release calls: the single production publisher. Main publishing runs independently of component checks, not after them. |
| `integration-test.yml` | Manual dispatch: publish a dedicated build, then run both complete real-data integrations. A reusable release call tests the supplied generated commit without republishing. No scheduled, push or PR integration runs. |
| `release-build.yml` | Manual dispatch with required `vX.Y.Z` version: publish first, then run component checks and both integrations in parallel against the exact generated commit. |

For manual runs, open the corresponding workflow under **Actions → Run workflow**
and select a **source branch**. Its selected immutable source commit is built.
Do not edit or dispatch from generated output branches:

| Build | Generated branch | Local production image tag / generated `_viash.yaml` version |
| --- | --- | --- |
| Main push | `main_build` | `main_build` |
| Standalone integration | `integration_build` | `integration_build` |
| Release | `release` | Supplied `vX.Y.Z`, e.g. `v0.1.0` |

Images use `ghcr.io/sunneyxie-lab/cpl-sc-foodie/<component>:<version>`, preserving
namespaced component paths. Biobox dependencies retain their upstream images.
Generated branches contain source and the complete production `target/`, not
fixtures or test-only targets/images. Integration checkers are built locally,
never pushed; production containers are pulled from published references.
Each integration downloads and validates common fixtures before its workflow
fixtures anonymously. See the existing [fixture and integration instructions](../resources_test_scripts/README.md)
for local usage and scientific coverage.

### Consume and tag a release

`main_build`, `integration_build` and `release` move on later builds. Version tags
are manual: no workflow creates Git tags or GitHub Release objects, or removes
historical release branches/tags. After a successful release run, inspect its
recorded generated `published_sha` and tag that exact commit:

```bash
git fetch origin release
git tag v0.1.0 <published_sha>
git push origin refs/tags/v0.1.0
```

Never tag source `main` or a later moving `release` tip instead. To consume it:

```bash
git fetch origin refs/tags/v0.1.0:refs/tags/v0.1.0
git checkout --detach v0.1.0
repo="$(pwd)"
```

Run `"$repo/target/nextflow/workflows/<name>/main.nf"` with the analysis arguments
in the [README](../README.md); do not rebuild the published tree. Keeping the tag preserves that
generated tree, not overwritten Docker image bytes. For stable consumption,
never reuse its release image version.

### Repository, GHCR and runner prerequisites

The public-runner design requires the user to make the repository public.
That disclosure includes existing code, branches and history; it is a user-owned
decision, not an automated visibility change or a sensitive-history audit.

Publishing uses the workflow's `GITHUB_TOKEN`, not a PAT or local token.
For each existing package under `ghcr.io/sunneyxie-lab/cpl-sc-foodie/`, grant
`sunneyxie-lab/cpl-sc-foodie` **Write** access under **Package settings → Manage
Actions access**. Allow public packages in the organization's package settings.
Set every production package's visibility to **Public**, including `fditools`,
`seqkit`, `trim_galore` and `restore_consensus_tags`.

New GHCR packages start private: a first run may push images but fail anonymous
manifest checks. Make those packages public, then **Re-run all jobs**.
The generated branch is published only after every image passes its anonymous
check. Public package visibility cannot be changed back to private. Workflows
never change repository or package visibility automatically.

Standard public Linux runners provide 4 CPUs, 16 GB RAM and 14 GB SSD; private
runners provide only 2 CPUs/8 GB and do not establish intended integration
acceptance. CI retains 4-CPU/13-GB process requests and `process.maxForks = 1`.
Disk, total memory and runtime sufficiency remain unmeasured. Each scientific
integration command has a 120-minute budget, not a measured runtime guarantee.

### Failures, retries and concurrency

Publication is not atomic and has no automatic rollback or retry:

- **Before publication:** invalid inputs, generation or component-discovery
  failures leave remote state unchanged.
- **During image pushes or anonymous checks:** some tags may already be replaced,
  while the generated branch stays unchanged. Cancellation, timeout or a lost
  runner can leave the same partial state. Inspect Actions logs and registry
  state, resolve access/capacity issues, then **Re-run all jobs**.
- **During branch publication:** images already exist. A failed push leaves the
  previous branch, but connection loss after server acceptance may leave the new
  branch despite a failed job. Inspect the remote commit before a full rerun.
  Failure to resolve the generated SHA or write the summary after a push does
  not undo publication.
- **After publication:** failed/cancelled component or integration checks,
  fixture validation failures and runner exhaustion leave branches and images
  published. Publication alone does not certify a release.

Publishing retries, including main, require **Re-run all jobs**; stale partial
retries are rejected before remote mutation. Standalone integration retries
also require a fresh full publish-and-test attempt, not failed-jobs-only retries.
Release test-only retries can reuse the recorded generated SHA and versioned
images. Full reruns can republish tags, including release versions: do not reuse
a version if its image bytes must stay stable.

Publishers serialize by mode; all release versions share one publisher lock.
Standalone integration holds one shared lock across publication and both tests,
regardless of source branch. Releases serialize publication and both check
suites. In-progress runs are not automatically cancelled, but GitHub's default
concurrency queue can replace older pending runs; not every queued run will run.
Serialization prevents workflow publishers overlapping, not external consumers
pulling during publication. Shared main/integration tags remain mutable.
Status and failure notifications use ordinary GitHub Actions statuses and logs.

