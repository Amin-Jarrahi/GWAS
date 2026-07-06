.PHONY: install install-dev run test lint clean help

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install:  ## Install package (editable)
	pip install -e .

install-dev:  ## Install with dev/test dependencies
	pip install -e ".[dev]"

run:  ## Run pipeline with default settings
	python scripts/run_pipeline.py

run-no-dowhy:  ## Run pipeline without DoWhy
	python scripts/run_pipeline.py --no_dowhy

test:  ## Run test suite
	pytest tests/ -v --tb=short

lint:  ## Lint source code with ruff
	ruff check src/ tests/ scripts/

format:  ## Auto-format with ruff
	ruff format src/ tests/ scripts/

typecheck:  ## Run mypy type checking
	mypy src/genomic_causal/

clean:  ## Remove caches and generated results
	rm -rf results/*.xlsx results/*.png
	rm -rf __pycache__ src/**/__pycache__ tests/**/__pycache__
	rm -rf .pytest_cache .mypy_cache .ruff_cache
	rm -rf *.egg-info dist build
