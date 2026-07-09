PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
ITENSOR_DIR := $(PROJECT_ROOT)/vendor/ITensor

.PHONY: help itensor check-itensor ising-builtin ising-manual ffd-manual ffd-production clean

help:
	@echo "Targets:"
	@echo "  make itensor        Build the vendored ITensor library"
	@echo "  make check-itensor  Verify the expected ITensor paths exist"
	@echo "  make ising-builtin  Build the library-assisted Ising tutorial"
	@echo "  make ising-manual   Build the manual Ising tutorial"
	@echo "  make clean          Remove local tutorial build artifacts"
	@echo "  make ffd-manual  	Build the FFD tutorial"
	@echo "  make ffd-production Build the production FFD TEBD runner"

check-itensor:
	@test -d "$(ITENSOR_DIR)" || (echo "Missing $(ITENSOR_DIR). Clone ITensor first." && false)
	@test -f "$(ITENSOR_DIR)/options.mk" || (echo "Missing $(ITENSOR_DIR)/options.mk. Copy options.mk.sample first." && false)

itensor: check-itensor
	$(MAKE) -C "$(ITENSOR_DIR)" -j1

ising-builtin: check-itensor
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ising_builtin"

ising-manual: check-itensor
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ising_manual"

ffd-manual:
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ffd_manual"

ffd-production:
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ffd_production"

clean:
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ising_builtin" clean
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ising_manual" clean
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ffd_manual" clean
	$(MAKE) -C "$(PROJECT_ROOT)/tutorials/ffd_production" clean
