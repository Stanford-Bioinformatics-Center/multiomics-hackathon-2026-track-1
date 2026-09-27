# Makefile for the STRING human network pipeline.
# Uses the project-local virtualenv in .venv.

PY := .venv/bin/python
PIP := .venv/bin/pip
CONFIG := config.yaml
THRESHOLD ?= 700

.PHONY: help setup data build build-500 compare test tune inspect clean \
        mnet-download mnet-metabolites mnet-lipids mnet-rhea mnet-assemble \
        mnet-qc mnet-test mnet-all mnet-ptm

# mnet is a src/ package; PYTHONPATH=src guarantees it is importable.
MNET_PY := PYTHONPATH=src $(PY)

help:
	@echo "Targets:"
	@echo "  make setup    - create .venv and install requirements"
	@echo "  make data     - download STRING + UniProt source files"
	@echo "  make build    - build the default output parquet (THRESHOLD?=700)"
	@echo "  make build-500- build the >= 500 output parquet"
	@echo "  make compare  - diff the build against the reference"
	@echo "  make tune     - grid-search build options against the reference"
	@echo "  make inspect  - profile the reference parquet"
	@echo "  make test     - run pytest"
	@echo "  make clean    - remove processed outputs"

setup:
	python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt
	$(PIP) install -e .

data:
	$(PY) -m string_network download --config $(CONFIG)

build:
	$(PY) -m string_network build --config $(CONFIG) --threshold $(THRESHOLD)

build-500:
	$(MAKE) build THRESHOLD=500

compare:
	$(PY) -m string_network compare --config $(CONFIG) --threshold $(THRESHOLD)

tune:
	$(PY) scripts/tune.py --config $(CONFIG)

inspect:
	$(PY) scripts/inspect_reference.py --config $(CONFIG)

test:
	$(PY) -m pytest -q

clean:
	rm -f data/processed/*.parquet

# --- Metabolite extension (mnet) ---------------------------------------------
mnet-download:
	$(MNET_PY) -m mnet download --config $(CONFIG)

mnet-metabolites:
	$(MNET_PY) -m mnet metabolites --config $(CONFIG)

mnet-lipids:
	$(MNET_PY) -m mnet lipids --config $(CONFIG)

mnet-rhea:
	$(MNET_PY) -m mnet rhea --config $(CONFIG)

mnet-assemble:
	$(MNET_PY) -m mnet assemble --config $(CONFIG)

mnet-qc:
	$(MNET_PY) -m mnet qc --config $(CONFIG)

mnet-test:
	$(MNET_PY) -m pytest tests/mnet -q

mnet-all: mnet-download mnet-metabolites mnet-lipids mnet-rhea mnet-assemble mnet-qc

mnet-ptm:
	$(MNET_PY) -m mnet ptm --config $(CONFIG)
