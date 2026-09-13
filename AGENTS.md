# Repository instructions

These instructions apply to all future work in this repository unless a task branch contains stricter task-specific instructions.

## Scientific task selection

Do not infer the active scientific task from stale historical branches or old prompts. The user will name the branch and prompt for each task. Read that task prompt first and obey its scientific scope and compute limits.

## GitHub publication workflow for ChatGPT Work

This repository is connected to GitHub through the authenticated GitHub integration. When operating in ChatGPT Work, use that connected GitHub integration for repository writes.

- Do not require shell `git push`, `gh auth`, SSH keys, credential helpers, PATs, or interactive Git authentication in order to complete or publish a task.
- Shell Git may be used locally for `status`, `diff`, inspection, and ordinary working-tree operations, but publication must not depend on shell credentials.
- Create or update branches, commit text/source/test outputs, and open or update pull requests through the connected GitHub integration when available.
- Before long computation, confirm that the requested remote branch exists and that the connected GitHub integration can read it.
- After computation, publish the completed source changes, focused tests, and compact UTF-8 scientific outputs to the requested branch through the connected GitHub integration.
- Do not use manual GitHub object construction, Base64 upload loops, archive upload, or hundreds of individual binary blobs as a workaround.
- Do not publish large binary trajectory collections unless the user explicitly requests them. Prefer compact CSV/JSON/Markdown summaries and source/test files.
- If the connected GitHub integration itself is unavailable, preserve local outputs and report the publication failure once. Do not rerun scientific calculations because publication failed.
- Never merge to `main` unless the user explicitly requests a merge.

## Repository discipline

Do not delete historical repository material merely to create a smaller working tree. If a reduced computation workspace is useful, keep it on a dedicated branch and do not merge that branch wholesale into `main`.

Do not edit `archive/` unless a task explicitly requires it. Preserve scientific provenance and distinguish experimental facts, model assumptions, calibration choices, and numerical diagnostics.
