##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic mesh tools project.
# Targets provided our detailed below...
#
# all: (default) Complete build and test the mesh tools.
#      Intended for interactive use.
# build: Build all executables and any supporting files
# test: Run test battery including unit tests and others
# clean: Delete all final products and working files
#
# The following variables may be specified to modify the build process...
#
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
# PROFILE: Set to a string representing a package of compiler options.
#          Potential profiles are 'full-debug', 'fast-debug' and 'production'.
#          Default: fast-debug
# LINK_TYPE: Either 'static' or 'dynamic'.
#            Default: dynamic
#
##############################################################################

PROJECT_NAME = mesh_tools
PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod planet_constants_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf hdf5 yaxt yaxt_c

PROJECT_DIR       := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR   ?= ..
export PATH       := $(PATH):$(abspath $(ROOT_DIR)/GPL/bin)
export PYTHONPATH := $(PYTHONPATH):$(abspath $(ROOT_DIR)/GPL/lib/python)

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build unit-tests
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk

export FFLAGS  += $(OPENMP_ARG)
export LDFLAGS += $(OPENMP_ARG)
ANY_UNIT_TESTS := $(shell if [ -d unit-test ] ; then \
                            find ./unit-test -name *test.pf ;\
                          fi)
ANY_INTEGRATION_TESTS := $(shell if [ -d integration-test ] ; then \
                                   find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \; ; \
                                 fi )

##############################################################################
# Launch gscan to monitor suites from this suite-group
#
.PHONY: launch-suite-gscan
launch-suite-gscan: gscan_processes := $(shell ps --no-headers -o command -C cylc-gscan)
launch-suite-gscan: ALWAYS
	$(Q)-if [[ "$(gscan_processes)" != *"--name=$(SUITE_GROUP_NAME)"* ]]; then \
          cylc gscan  $(DOUBLE_VERBOSE_ARG) --name=$(SUITE_GROUP_NAME) &           \
          usleep 1                                                                ;\
        fi

##############################################################################

.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Build
#
.PHONY: build
ifeq "$(PROFILE)" "full-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_INIT) \
            $(FFLAGS_RUNTIME) $(FFLAGS_NO_OPTIMISATION) $(STANDARDS_CHECK_ARG)
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) \
            $(FFLAGS_SAFE_OPTIMISATION) $(STANDARDS_CHECK_ARG)
else ifeq "$(PROFILE)" "production"
build: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_WARNINGS) $(FFLAGS_RISKY_OPTIMISATION)
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: PROJECT      = mesh_tools
build: WORKING_DIR := $(WORKING_DIR)/mesh_tools
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: ALWAYS
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=$(PROJECT) \
	          SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	          BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"


##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export PROJECT      = mesh_tools
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-test
unit-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
unit-tests: export PROGRAMS     = mesh_tools_unit_tests
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS      += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
ifeq ($(strip $(ANY_UNIT_TESTS)),)
unit-tests:
	$(call MESSAGE,No unit-tests to process.)
else
unit-tests:run-unit-tests
	$(Q)echo > /dev/null
endif

.PHONY: generate-unit-tests
generate-unit-tests: ALWAYS
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=$(PROJECT) \
	          SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk \
                  PROJECT=mesh_tools \
                  SOURCE_DIR=unit-test \
                  WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"mesh tools work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"mesh tools binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(call MESSAGE,Removing,"mesh tools tests")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
