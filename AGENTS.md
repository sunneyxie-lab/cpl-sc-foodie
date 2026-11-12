# AGENTS.md

Viash pipeline project based on [openpipeline_incubator](https://github.com/openpipelines-bio/openpipeline_incubator),
the OpenPipeline ecosystem's sample project. Ecosystem docs (draft):
[openpipeline_docs](https://github.com/openpipelines-bio/openpipeline_docs).
Components wrap CLI tools; workflows are Nextflow pipelines composed from components.
Both are declared via `config.vsh.yaml` and built/tested with the `viash` CLI.

## Structure

```
_viash.yaml                          # viash project config (repositories, ...)
nextflow.config                      # root Nextflow config
resources_test_scripts/              # bash scripts that fetch test resources
src/
  _authors/                          # shared author YAML, __merge__'d into configs
  base/
    requirements/                    # shared requirement snippets
    <name>_argument.yaml             # shared argument, e.g. h5_compression_argument.yaml
  _utils/                            # shared scripts, e.g. test_helpers.sh
  components/
    <tool>/                          # single wrapper; no `namespace:`
    <tool>/<tool>_<subcommand>/      # one per subcommand; `namespace: <tool>`
    dataflow/<name>/                 # split, merge, concatenate, ...; `namespace: dataflow`
  workflows/
    utils/                           # shared Nextflow configs, e.g. labels.config
    <pipeline_name>/
      config.vsh.yaml
      integration_test.sh
      main.nf
      nextflow.config
      test.nf
    test_workflows/<pipeline_name>_test/   # optional `scope: test` component checking the pipeline's output
```

Start a tool top-level; switch to the `<tool>_<subcommand>` layout only when its second
subcommand arrives. Group workflows into categories only once there are several pipelines.
`test_workflows/` stays in `src/workflows/` even then (openpipeline_incubator convention); only
its inner nesting mirrors the pipeline layout.

Avoid writing `.nf` files in `resources_test_scripts/`.

## Components

Reuse [biobox](https://github.com/viash-hub/biobox) components when available: register
biobox under `repositories:` in `_viash.yaml` and reference them from workflow configs
via `dependencies:` (`name: <namespace>/<name>`, `repository: biobox`).

Otherwise write a local component following biobox
[COMPONENT_DEVELOPMENT.md](https://github.com/viash-hub/biobox/blob/main/docs/COMPONENT_DEVELOPMENT.md) and the
docs it links.

For writing Python-based component scripts, follow the biobox
[Python Script Development Guide](https://github.com/viash-hub/biobox/blob/main/docs/PYTHON_SCRIPT_DEVELOPMENT.md):
`par` / `meta` dicts, subprocess handling, and when to pick Python over bash.

For setting up Docker engines and managing dependencies in components, follow the biobox
[Docker Guide](https://github.com/viash-hub/biobox/blob/main/docs/DOCKER_GUIDE.md).

For one-to-many/many-to-one steps (per-sample vs. all-samples tools) or per-sample
parallelism, add dataflow components.

## Workflows

Follow openpipeline's
[Create a pipeline](https://github.com/openpipelines-bio/openpipeline_docs/blob/main/contributing/creating-workflows.qmd)
guide, except with no category: `src/workflows/<pipeline_name>/`, `namespace: workflows`.

## Testing

Use synthetic data for component tests and real data for integration tests.

Do not require permanent companion tests for fixture-preparation scripts under `resources_test_scripts/`; use smoke checks and workflow integration unless dedicated regression coverage is justified.

Components: follow biobox
[Test conventions](https://github.com/viash-hub/biobox/blob/main/docs/COMPONENT_DEVELOPMENT.md#testing)
and its commands.

Pipeline integration test (full Nextflow pipeline), from the repo root:

```bash
src/workflows/<pipeline_name>/integration_test.sh
```

If needed, write a `test_workflows/` component; `test.nf` pipes the workflow's output into it.

## Docs

- [Agent setup](.agents/docs/onboarding.md) — read when wiring up an agent or changing the agent setup.

## Plans

Planning docs (plans, requirements, investigations) go in `.agents/plans/` (gitignored).
