# Publishing Guide (English)

This guide covers everything needed to build **opticspy-research**, upload it to PyPI, and
set up **automatic publishing** so that every new GitHub Release is pushed to PyPI by itself
-- no manual `twine upload` needed after the first setup.

The Python package importable as `import opticspy` is distributed on PyPI under the name
**`opticspy-research`** (the plain name `opticspy` was already taken by an unrelated project).

---

## 0. One-time prerequisites

1. A [PyPI](https://pypi.org/account/register/) account (and, recommended, a
   [TestPyPI](https://test.pypi.org/account/register/) account for dry runs).
2. A GitHub repository containing this project (push the contents of this zip to a new repo).
3. Before your first real upload, **edit `pyproject.toml`**: replace every `OWNER` in the
   `[project.urls]` section with your actual GitHub username/organization, and update the
   `authors` field if you want a named maintainer instead of "OpticsPy Contributors".
4. Python 3.9+ and `pip install build twine` locally (already included via `pip install -e .[dev]`).

---

## 1. Build and check the package locally

```bash
cd opticspy-research/          # the folder containing pyproject.toml
python -m pip install --upgrade build twine
python -m build                # creates dist/*.whl and dist/*.tar.gz
twine check dist/*             # validates metadata/long-description rendering
```

Run the test suite first to make sure everything passes:
```bash
pip install -e .[dev]
pytest tests/ -q
```

---

## 2. First upload: manual, via TestPyPI then PyPI

It's good practice to upload to TestPyPI first to make sure everything renders and installs
correctly before publishing for real.

```bash
# 1) Upload to TestPyPI
twine upload --repository testpypi dist/*
# (enter your TestPyPI API token when prompted, username: __token__)

# 2) Install from TestPyPI in a clean venv to verify
python -m venv /tmp/test_env && source /tmp/test_env/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ opticspy-research
python -c "import opticspy; print(opticspy.__version__)"
deactivate

# 3) If everything looks right, upload to the real PyPI
twine upload dist/*
# (username: __token__, password: your PyPI API token)
```

Getting an API token: PyPI account settings -> "API tokens" -> "Add API token". Scope it to
the `opticspy-research` project (after the very first upload, which needs an account-wide
token; you can narrow the scope afterwards).

After this first manual upload, `pip install opticspy-research` works for everyone.

---

## 3. Automatic publishing on every future release (GitHub Actions)

This repository already includes `.github/workflows/publish.yml`. It runs whenever you
publish a **GitHub Release**, builds the package, and uploads it to PyPI automatically using
**Trusted Publishing** (OpenID Connect) -- the current PyPI-recommended method, which needs
**no secret API token stored in GitHub** at all.

### 3.1 One-time setup: link GitHub to PyPI (Trusted Publishing)

1. Go to your project on PyPI: `https://pypi.org/manage/project/opticspy-research/settings/publishing/`
   (you must have done the first manual upload from Section 2 so the project exists).
2. Under "Trusted Publishers", click "Add a new publisher" and fill in:
   - Owner: your GitHub username/organization
   - Repository name: your repo name (e.g. `opticspy-research`)
   - Workflow name: `publish.yml`
   - Environment name: `pypi`
3. Save. That's it -- no tokens to copy anywhere.

### 3.2 How a release triggers a publish

1. Bump the version number in **`opticspy/__init__.py`** (`__version__ = "2.2.0"`, etc.) --
   this is the single source of truth; `pyproject.toml` reads it automatically.
2. Add a matching entry at the top of **`CHANGELOG.md`**.
3. Commit, then tag and push:
   ```bash
   git add -A && git commit -m "Release v2.2.0"
   git tag v2.2.0
   git push && git push --tags
   ```
4. On GitHub, go to "Releases" -> "Draft a new release", pick the tag you just pushed
   (e.g. `v2.2.0`), write release notes (or copy from CHANGELOG.md), and click
   "Publish release".
5. Publishing the release automatically triggers `.github/workflows/publish.yml`, which
   builds the package and uploads it to PyPI within a minute or two. No further action needed.

### 3.3 Fallback: API-token method (if you prefer not to use Trusted Publishing)

If you'd rather use a classic token instead of OIDC trusted publishing:
1. Create a PyPI API token scoped to the project.
2. In your GitHub repo: Settings -> Secrets and variables -> Actions -> "New repository
   secret", name it `PYPI_API_TOKEN`, paste the token.
3. Replace the `publish` job's steps in `.github/workflows/publish.yml` with:
   ```yaml
   - uses: pypa/gh-action-pypi-publish@release/v1
     with:
       password: ${{ secrets.PYPI_API_TOKEN }}
   ```
   (and drop the `permissions: id-token: write` block, which is only needed for OIDC).

---

## 4. Continuous integration (tests on every push/PR)

`.github/workflows/tests.yml` runs `pytest` across Python 3.9-3.12 on every push and pull
request to `main`, so regressions are caught before a release is even tagged.

---

## 5. Quick checklist for every new release

- [ ] Code changes merged to `main`, all tests green.
- [ ] `opticspy/__init__.py`: `__version__` bumped.
- [ ] `CHANGELOG.md`: new entry added.
- [ ] `git tag vX.Y.Z && git push --tags`
- [ ] GitHub "Draft a new release" from that tag -> "Publish release".
- [ ] Watch the Actions tab for the `publish.yml` run to go green.
- [ ] `pip install --upgrade opticspy-research` and confirm the new version.
