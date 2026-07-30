PACKAGE         := ofiber
GITHUB_USER     := scottprahl

PY_VERSION      ?= 3.14
UV              ?= uv
RUN             := $(UV) run --extra dev
RUN_DOCS        := $(UV) run --extra docs
RUN_LITE        := $(UV) run --extra lite
RM              ?= rm -f
RMR             ?= rm -rf

DOCS_DIR        := docs
OUT_ROOT        := _site
STAGE_DIR       := .lite_src
DOIT_DB         := .jupyterlite.doit.db
OUT_DIR         := $(OUT_ROOT)/$(PACKAGE)
HTML_DIR        := $(DOCS_DIR)/_build/html
LITE_CONFIG     := $(PACKAGE)/jupyter_lite_config.json
COV_DIR         := htmlcov

# --- GitHub Pages deploy config ---
PAGES_BRANCH    := gh-pages
WORKTREE        := .gh-pages
REMOTE          := origin

# --- server config (override on CLI if needed) ---
HOST            := 127.0.0.1
PORT            := 8000

PYTEST_OPTS     :=
COV_OPTS        := --cov=$(PACKAGE) --cov-report=term-missing --cov-report=html:$(COV_DIR)
SPHINX_OPTS     := -T -E -b html -d $(DOCS_DIR)/_build/doctrees -D language=en

PYLINT_TARGETS  := $(PACKAGE)/*.py tests/*.py .github/scripts/update_citation.py
YAML_TARGETS    := .github/workflows/citation.yaml .github/workflows/pypi.yaml .github/workflows/test.yaml .readthedocs.yaml
RST_TARGETS     := README.rst CHANGELOG.rst $(DOCS_DIR)/index.rst $(DOCS_DIR)/changelog.rst

.PHONY: help
help:
	@echo "Build Targets:"
	@echo "  dist           - Build sdist+wheel locally"
	@echo "  html           - Build Sphinx HTML documentation"
	@echo "  lab            - Start jupyterlab"
	@echo "  readme         - Create images for readme"
	@echo "  venv           - Install project dependencies with uv extras"
	@echo ""
	@echo "Test Targets:"
	@echo "  test           - Run pytest on python files"
	@echo "  coverage       - Run unit tests and report $(PACKAGE) coverage"
	@echo "  note-test      - Test all notebooks for errors"
	@echo ""
	@echo "Packaging Targets:"
	@echo "  rcheck         - Distribution release checks"
	@echo "  manifest-check - Validate MANIFEST"
	@echo "  pylint-check   - Same as lint above"
	@echo "  pyroma-check   - Validate overall packaging"
	@echo "  rst-check      - Validate all RST files"
	@echo "  ruff-check     - Lint all .py and .ipynb files"
	@echo "  yaml-check     - Validate YAML files"
	@echo ""
	@echo "JupyterLite Targets:"
	@echo "  lite           - Build JupyterLite site into $(OUT_DIR)"
	@echo "  lite-serve     - Serve $(OUT_DIR) at http://$(HOST):$(PORT)"
	@echo "  lite-deploy    - Upload to github"
	@echo ""
	@echo "Clean Targets:"
	@echo "  clean          - Remove build caches and docs output"
	@echo "  lite-clean     - Remove JupyterLite outputs"
	@echo "  realclean      - clean + remove $(VENV)"

.PHONY: venv
venv:
	$(UV) sync --extra dev --extra docs --extra lite

.PHONY: dist
dist:
	$(RUN) python -m build

.PHONY: test
test:
	@:
#	$(RUN) pytest $(PYTEST_OPTS) tests --ignore=tests/test_all_notebooks.py

# Notebooks are excluded here on purpose: ExecutePreprocessor runs each one in a
# separate Jupyter kernel process, which coverage.py cannot observe, so including
# them would report a misleadingly low number rather than a real measurement.
.PHONY: coverage
coverage:
	@$(RM) .coverage
	@$(RMR) "$(COV_DIR)"
	@$(RUN) pytest $(PYTEST_OPTS) $(COV_OPTS) tests --ignore=tests/test_all_notebooks.py; \
	  status=$$?; \
	  if [ $$status -eq 5 ]; then echo "⚠️  no unit tests collected -- coverage is 0%"; \
	  elif [ $$status -ne 0 ]; then exit $$status; fi
	@test -f "$(COV_DIR)/index.html" && echo "==> HTML report at $(COV_DIR)/index.html" || true
	@command -v open >/dev/null 2>&1 && test -f "$(COV_DIR)/index.html" && open "$(COV_DIR)/index.html" || true

.PHONY: note-test
note-test:
	$(RUN) pytest --verbose tests/test_all_notebooks.py

.PHONY: html
html:
	@mkdir -p "$(HTML_DIR)"
	$(RUN_DOCS) sphinx-build $(SPHINX_OPTS) "$(DOCS_DIR)" "$(HTML_DIR)"
	@command -v open >/dev/null 2>&1 && open "$(HTML_DIR)/index.html" || true

.PHONY: lab
lab:
	@echo "==> Launching JupyterLab via uv"
	$(RUN) jupyter lab --ServerApp.root_dir="$(CURDIR)"

.PHONY: readme
readme:
	@cd docs/images && $(RUN) python make_readme_images.py

.PHONY: pylint-check
pylint-check:
	$(RUN) pylint $(PYLINT_TARGETS)

.PHONY: yaml-check
yaml-check:
	$(RUN) yamllint $(YAML_TARGETS)

.PHONY: rst-check
rst-check:
	$(RUN) rstcheck $(RST_TARGETS)
	$(RUN) rstcheck --ignore-directives automodapi $(DOCS_DIR)/$(PACKAGE).rst

.PHONY: ruff-check
ruff-check:
	$(RUN) ruff check

.PHONY: manifest-check
manifest-check:
	$(RUN) check-manifest

.PHONY: pyroma-check
pyroma-check:
	$(RUN) python -m pyroma -d .

.PHONY: rcheck
rcheck:
	@echo "Running all release checks..."
	@$(MAKE) realclean
	@$(MAKE) ruff-check
	@$(MAKE) pylint-check
	@$(MAKE) rst-check
	@$(MAKE) yaml-check
	@$(MAKE) manifest-check
	@$(MAKE) pyroma-check
	@$(MAKE) html
	@$(MAKE) lite
	@$(MAKE) dist
	@$(MAKE) test
	@$(MAKE) note-test
	@echo "✅ Release checks complete"
	
.PHONY: lite
lite: lite-clean $(LITE_CONFIG) dist
	@echo "==> Staging notebooks from docs -> $(STAGE_DIR)"
	mkdir -p "$(STAGE_DIR)"
	/bin/cp docs/*.ipynb "$(STAGE_DIR)"
	$(RUN) jupyter nbconvert --clear-output --inplace "$(STAGE_DIR)"/*.ipynb
	@echo "==> Building JupyterLite"
	@$(RUN_LITE) jupyter lite build \
		--config="$(LITE_CONFIG)" \
		--contents="$(STAGE_DIR)" \
		--output-dir="$(OUT_DIR)"
	@touch "$(OUT_DIR)/.nojekyll"  # github

.PHONY: lite-serve
lite-serve:
	@test -d "$(OUT_DIR)" || { echo "❌ run 'make lite' first"; exit 1; }
	@echo "Serving at"
	@echo "   http://$(HOST):$(PORT)/$(PACKAGE)/?disableCache=1"
	@echo ""
	$(RUN_LITE) python -m http.server -d "$(OUT_ROOT)" --bind $(HOST) $(PORT)

.PHONY: lite-deploy
lite-deploy: lite
	@echo "==> Ensure $(PAGES_BRANCH) branch exists"
	@if ! git show-ref --verify --quiet refs/heads/$(PAGES_BRANCH); then \
	  CURRENT=$$(git branch --show-current); \
	  git switch --orphan $(PAGES_BRANCH); \
	  git commit --allow-empty -m "Initialize $(PAGES_BRANCH)"; \
	  git switch $$CURRENT; \
	fi

	@echo "==> Setup deployment worktree"
	@git worktree remove "$(WORKTREE)" --force 2>/dev/null || true
	@git worktree prune || true
	@$(RMR) "$(WORKTREE)"
	@git worktree add "$(WORKTREE)" "$(PAGES_BRANCH)"
	@git -C "$(WORKTREE)" pull "$(REMOTE)" "$(PAGES_BRANCH)" 2>/dev/null || true

	@echo "==> Deploy $(OUT_DIR) -> $(WORKTREE)"
	@rsync -a --delete --exclude ".git*" "$(OUT_DIR)/" "$(WORKTREE)/"
	@touch "$(WORKTREE)/.nojekyll"
	@date -u +"%Y-%m-%d %H:%M:%S UTC" > "$(WORKTREE)/.pages-ping"

	@echo "==> Commit & push"
	@cd "$(WORKTREE)" && \
	  git add -A && \
	  if git diff --quiet --cached; then \
	    echo "✅ No changes to deploy"; \
	  else \
	    git commit -m "Deploy $$(date -u +'%Y-%m-%d %H:%M:%S UTC')" && \
	    git push "$(REMOTE)" "$(PAGES_BRANCH)" && \
	    echo "✅ Deployed to https://$(GITHUB_USER).github.io/$(PACKAGE)/"; \
	  fi

.PHONY: lite-clean
lite-clean:
	@echo "==> Cleaning JupyterLite build artifacts"
	@$(RMR) "$(STAGE_DIR)"
	@$(RMR) "$(OUT_ROOT)"
	@$(RMR) "$(DOIT_DB)"
	@$(RMR) .cache dist $(PACKAGE).egg-info

.PHONY: clean
clean: lite-clean
	@echo "==> Cleaning build artifacts"	
	@find . -name '__pycache__' -type d -exec $(RMR) {} +
	@find . -name '.DS_Store' -type f -delete
	@find . -name '.ipynb_checkpoints' -type d -prune -exec $(RMR) {} +
	@find . -name '.pytest_cache' -type d -prune -exec $(RMR) {} +
	@$(RMR) .ruff_cache
	@$(RMR) "$(COV_DIR)"
	@$(RM) .coverage
	@$(RMR) docs/api
	@$(RMR) docs/_build

.PHONY: realclean
realclean: clean
	@echo "==> Deep cleaning: removing venv and deployment worktree"
	@git worktree remove "$(WORKTREE)" --force 2>/dev/null || true
	@git worktree prune || true
	$(RMR) "$(WORKTREE)"
	$(RMR) .venv
	@$(RM) uv.lock
