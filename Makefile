.PHONY: test

test:
	python3 -m pytest -v --flake8 --pylint tests/MetaPhinder_test.py MetaPhinder.py
