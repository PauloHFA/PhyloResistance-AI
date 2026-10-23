.PHONY: install test lint pipeline api clean

install:
	uv sync --all-extras

test:
	uv run pytest -v

lint:
	uv run ruff check src tests

pipeline:
	uv run python -m phyloresistance.cli run-all --output-dir data/processed --reports-dir reports

api:
	uv run uvicorn phyloresistance.api.app:app --host 0.0.0.0 --port 8000 --reload

clean:
	rm -rf .pytest_cache .ruff_cache .mypy_cache tmp_uploads
