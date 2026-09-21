# iEEGLAB: notes for Claude

Read these first:
- `.claude/LESSONS.md`: hard-won technical lessons (newest last).
- `.claude/HANDOFF.md`: state of the work as of the last session, what is pending, where data live.
- `TODO.md` and `CHANGELOG.md` ("Unreleased" section) for what changed and what is open.
- Latest analysis write-up: `worklog/2026-09-20-native-rate-numbers.md` (read section 8,
  corrections, before quoting any number from sections 1-7).

Conventions:
- Tests are compute-only (`runtests('tests')`); MATLAB `-batch` hangs on figures on this machine.
- Prose written for Cedric or his collaborators: no em dashes; APA references; concise.
- `fix/carla-and-open-issues` (including the 2026-09-19..21 Cowork work) was merged into `main`
  on 2026-09-21. Review `git status` before committing.
