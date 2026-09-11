# Contributing to simulation_coarse

Thank you for your interest in contributing!

## How to contribute

1. Fork the repository and create a feature branch (`git checkout -b feature/my-feature`).
2. Make your changes with clear, focused commits.
3. Ensure all existing tests pass and add tests for new functionality.
4. Open a pull request against `main` with a concise description of the change.
5. A maintainer will review your PR; please address feedback promptly.

We welcome bug reports, feature requests, and documentation improvements via GitHub Issues.

## Dev setup

No build step is required — the project is plain JavaScript (ES modules) plus a static `index.html`.

```bash
# clone and serve locally
git clone https://github.com/<org>/simulation_coarse.git
cd simulation_coarse
python3 -m http.server
# open http://localhost:8000
```

Prerequisites: a modern browser, Node.js 18+ (for tests and syntax checks), and Python 3 (for the dev server).

## Tests

```bash
# syntax check
node --check src/*.js src/physics/*.js

# full test suite
npm test
# or directly:
node tests/test_all.js
node tests/test_gb_fd.js
node tests/test_forces_fd.js
```

CI (`.github/workflows/check.yml`) runs the same commands on every push and pull request.

## Good first issue

New contributors — look for issues labeled `good first issue` on the GitHub Issues page. These are scoped, well-described tasks ideal for onboarding:

- Documentation fixes and typo corrections
- Adding or improving unit tests
- Small UI/UX polish in `src/ui.js` or `src/viewer.js`
- Extending `data/manifest.json` entries for new reference structures

If you would like to be assigned, comment on the issue and a maintainer will triage it. Feel free to open a new issue and label it `good first issue` if you spot a suitable task.
