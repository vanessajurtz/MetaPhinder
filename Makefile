.PHONY: test

test:
	python3 -m pytest -v --flake8 --pylint --mypy MetaPhinder.py tests/MetaPhinder_test.py
