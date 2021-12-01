.PHONY: test

test:
	python3 -m pytest -v --flake8 --pylint MetaPhinder.py tests/MetaPhinder_test.py
