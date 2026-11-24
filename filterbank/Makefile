CARGO ?= cargo
NVCC ?= $(shell command -v nvcc 2>/dev/null)
PY ?= python3

.PHONY: all rust cuda test clean

all: rust cuda

rust:
	$(CARGO) build --release

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

clean:
	$(CARGO) clean
	-$(MAKE) -C cuda clean
