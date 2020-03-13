##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_NAME = catalyst_demo
export PROJECT_DIR_NAME = catalyst_demo

PROFILE ?= fast-debug

export IGNORE_DEPENDENCIES = netcdf MPI yaxt pfunit_mod xios mod_wait
export EXTERNAL_DYNAMIC_LIBRARIES = yaxt yaxt_c netcdff netcdf hdf5 \
                                    catalyst_adaptor vtkPVPythonCatalyst-pv5.5 \
                                    $(CXX_RUNTIME_LIBRARY)
export EXTERNAL_STATIC_LIBRARIES = xios

LFRIC_TARGET_PLATFORM ?= generic-x86
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR ?= ../..
export PATH        := $(PATH):$(abspath $(ROOT_DIR)/GPL/bin)
export PYTHONPATH  := $(PYTHONPATH):$(abspath $(ROOT_DIR)/GPL/lib/python)

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build unit-tests integration-tests
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-latex document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include $(ROOT_DIR)/infrastructure/build/fortran/$(FORTRAN_COMPILER).mk

export FFLAGS  := $(FFLAGS)  $(OPENMP_ARG)
export LDFLAGS := $(LDFLAGS) $(OPENMP_ARG)
ANY_UNIT_TESTS := $(shell if [ -d unit-test ] ; then \
                            find unit-test -name *test.pf ; \
                          fi )
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
# Documentation
#
.PHONY: document-uml
document-uml: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/uml
document-uml: SOURCE_DIR    = documentation/uml
document-uml: WORKING_DIR  := $(WORKING_DIR)/uml
document-uml: uml-documentation
	$(Q)echo > /dev/null

.PHONY: document-latex
document-latex: DOCUMENT_DIR  ?= $(PROJECT_DIR)/documents
document-latex: SOURCE_DIR     = documentation
document-latex: DOCUMENTS      = $(shell find $(SOURCE_DIR) -name '*.latex' -print)
document-latex: WORKING_DIR   := $(WORKING_DIR)/latex
document-latex: TEX_STUFF      = documentation/tex
document-latex: COMMON_FIGURES = documentation/common-figures
document-latex: ALWAYS
	$(Q)$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/latex.mk       \
                                 SOURCE_DIR=$(SOURCE_DIR)         \
	                         WORKING_DIR=$(WORKING_DIR)       \
	                         DOCUMENT_DIR=$(DOCUMENT_DIR)     \
	                         TEX_STUFF=$(TEX_STUFF)           \
	                         COMMON_FIGURES=$(COMMON_FIGURES) \
	                         DOCUMENTS="$(DOCUMENTS)"

.PHONY: document-api
document-api: PROJECT       = $(PROJECT_NAME)
document-api: DOCUMENT_DIR ?= $(PROJECT_DIR)/documents/api
document-api: CONFIG_DIR    = documentation
document-api: SOURCE_DIR    = source
document-api: WORKING_DIR  := $(WORKING_DIR)/api
document-api: api-documentation
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
build: BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/lfric-catalyst_demo/HEAD
build: PROGRAMS    := $(basename $(notdir $(shell find $(SOURCE_DIR) -maxdepth 1 -name '*.[Ff]90' -print)))
build: PROJECT      = $(PROJECT_NAME)
build: ALWAYS
	$(call MESSAGE,===============================)
	$(call MESSAGE,Extracting LFRic infrastructure)
	$(call MESSAGE,===============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk  \
	          SOURCE_DIR=$(ROOT_DIR)/infrastructure/source \
	          WORKING_DIR=$(WORKING_DIR)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/infrastructure/source \
	          WORKING_DIR=$(WORKING_DIR) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	          WORKING_DIR=$(WORKING_DIR) \
	          WITHOUT_PROGRAMS=1
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=gungho \
	          SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	          WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$(call MESSAGE,====================================)
	$(call MESSAGE,Psycloning Gungho dynamical core)
	$(call MESSAGE,====================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	          WORKING_DIR=$(WORKING_DIR) \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,==============================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,==============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR)
	$(call MESSAGE,==============================)
	$(call MESSAGE,Psycloning $(PROJECT_NAME))
	$(call MESSAGE,==============================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone.mk \
	          SOURCE_DIR=source \
		  WORKING_DIR=$(WORKING_DIR) \
                  OPTIMISATION_PATH=$(OPTIMISATION_PATH)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Processing namelist loaders for $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  PROJECT=$(PROJECT) \
	          SOURCE_DIR=source \
	          WORKING_DIR=$(WORKING_DIR) \
                  META_FILE_DIR=$(META_FILE_DIR) \
                  FIX_ENUMS_OPT=$(FIX_ENUMS_OPT)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$(Q)$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	         BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"


##############################################################################
# Unit tests
#
.PHONY: unit-tests
unit-tests: export ADDITIONAL_EXTRACTION = "$(ROOT_DIR)/infrastructure/source \
                                            $(ROOT_DIR)/gungho/source"
unit-tests: export BIN_DIR = $(PROJECT_DIR)/test
unit-tests: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests: export FFLAGS += $(FFLAGS_DEBUG) $(FFLAGS_NO_OPTIMISATION) $(FFLAGS_INIT)
unit-tests: export IGNORE_DEPENDENCIES += pfunit_mod
unit-tests: export META_FILE_DIR = rose-meta/lfric-catalyst_demo/HEAD
unit-tests: export PROGRAMS = $(PROJECT_NAME)_unit_tests
unit-tests: export PROJECT = $(PROJECT_NAME)
unit-tests: export SOURCE_DIR = source
unit-tests: export TEST_DIR = unit-test
unit-tests: export WORKING_DIR := $(WORKING_DIR)/unit-test
unit-tests: do-unit-tests
	$(Q)echo > /dev/null


##############################################################################
# Integration tests
#
.PHONY: integration-tests
integration-tests: export BIN_DIR      = $(PROJECT_DIR)/test
integration-tests: export PROJECT      = $(PROJECT_NAME)
integration-tests: export SOURCE_DIR   = source
integration-tests: export TEST_DIR     = integration-test
integration-tests: export WORKING_DIR := $(WORKING_DIR)/integration-test
integration-tests: do-integration-tests
	$(Q)echo > /dev/null


##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"$(PROJECT_NAME) work space")
	$(Q)-if [ $(WORKING_DIR) != *[\*]* ] && [ -d $(WORKING_DIR) ] ; then rm -r $(WORKING_DIR) ; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) documents")
	$(Q)-if [ -d documents ] ; then rm -r documents; fi
	$(call MESSAGE,Removing,"$(PROJECT_NAME) binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi
	$(Q)-if [ -d test ] ; then rm -r test ; fi
