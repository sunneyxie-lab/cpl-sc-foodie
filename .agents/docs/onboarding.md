# Agent setup

`AGENTS.md` is the canonical, agent-neutral entry file: a lean index linking to the docs in
`.agents/docs/`, loaded on demand. Other layout:

- `.agents/skills/<name>/` — shared project skills; each agent reaches them via a relative symlink
  in its own skills dir (unless it reads `.agents/skills` natively).
- `.agents/plans/` — planning docs; gitignored, so personal.

## Wiring agents

- **Claude Code** reads `CLAUDE.md`, not `AGENTS.md`: `CLAUDE.md` holds `@AGENTS.md` (plus any
  Claude-only lines below it); skills are linked as `.claude/skills/<name>` →
  `../../.agents/skills/<name>`.
- **Any other agent**: point its entry-file mechanism at `AGENTS.md` (an include, or a copy
  kept in sync where includes are unsupported), or re-run `/agentify-project` to research and
  wire it.

## Windows

Skill symlinks need Developer Mode (or admin) and `git config core.symlinks true` before
cloning; otherwise they check out as plain text files.
