CARGO ?= cargo
NVCC ?= $(shell command -v nvcc 2>/dev/null)
PY ?= python3

.PHONY: all rust dada cuda test clean

all: rust dada cuda

rust:
	$(CARGO) build --release

dada:
	$(MAKE) -C dada

cuda:
ifneq ($(NVCC),)
	$(MAKE) -C cuda
else
	@echo "nvcc not found: CPU backend only"
endif

test: all
	$(CARGO) test --release -q
	PYTHONPATH=python $(PY) tests/test_python_api.py
	$(PY) tests/test_sp_cli.py
	$(PY) tests/test_psrfits_convert.py
	$(PY) tests/test_periodic_cli.py
	$(PY) tests/test_toa_cli.py
	$(PY) tests/test_classify_realtime.py
	$(PY) tests/test_ring.py

clean:
	$(CARGO) clean
	-$(MAKE) -C dada clean
	-$(MAKE) -C cuda clean
