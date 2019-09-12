##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic infrastructure project.
# Targets provided our detailed below...
#
# all: (default) Complete build and test the infrastructure.
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

PROJECT_NAME = infrastructure

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod planet_constants_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = netcdff netcdf hdf5 yaxt yaxt_c

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..
export PATH        := $(PATH):$(abspath $(ROOT_DIR)/GPL/bin)
export PYTHONPATH  := $(PYTHONPATH):$(abspath $(ROOT_DIR)/GPL/lib/python)

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-api
	$(Q)echo > /dev/null

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

include build/lfric.mk

export FFLAGS  := $(FFLAGS) $(OPENMP_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)
ANY_UNIT_TESTS := $(shell if [ -d unit-test ] ; then \
                            find unit-test -name *test.pf ; \
                          fi )
ANY_INTEGRATION_TESTS := $(shell if [ -d integration-test ] ; then \
                                   find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \; ; \
                                 fi )


##############################################################################
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null


.PHONY: document-api
document-api: PROJECT       = infrastructure
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
	$(Q)echo > /dev/null


##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export PROJECT      = infrastructure
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: export BIN_DIR     ?= $(PROJECT_DIR)/test
unit-tests: export PROGRAMS     = infrastructure_unit_tests
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
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=unit-test \
                  WORKING_DIR=$(WORKING_DIR)
	$Qcp -r unit-test/data $(WORKING_DIR)/
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/pfunit.mk \
                  PROJECT=infrastructure \
                  SOURCE_DIR=unit-test \
                  WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export PROJECT        = infrastructure
integration-tests: export WORKING_DIR   := $(WORKING_DIR)/integration-tests
integration-tests: export BIN_DIR       ?= $(PROJECT_DIR)/test
ifeq ($(strip $(ANY_INTEGRATION_TESTS)),)
integration-tests:
	$(call MESSAGE,No integration-tests to process.)
else
integration-tests: $(addprefix integration-test-run/,$(basename $(shell find integration-test -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
	$(Q)echo > /dev/null
endif

.PHONY: generate-integration-tests
generate-integration-tests: ALWAYS
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=infrastructure \
                  SOURCE_DIR=source \
                  WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=integration-test/configuration_test \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=integration-test \
                  WORKING_DIR=$(WORKING_DIR)


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"infrastructure work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"infrastructure documents")
	$(Q)if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"infrastructure binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
