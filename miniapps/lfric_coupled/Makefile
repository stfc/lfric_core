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

export PROJECT_NAME = lfric_coupled
export PROJECT_DIR_NAME = lfric_coupled
PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ../..

export IGNORE_DEPENDENCIES += mod_oasis
export EXTERNAL_STATIC_LIBRARIES += psmile.MPI1 mct mpeu scrip

export INTERNAL_DEPENDENCIES = $(ROOT_DIR)/infrastructure \
                               $(ROOT_DIR)/components/driver \
                               $(ROOT_DIR)/components/science \
                               $(ROOT_DIR)/components/inventory \
                               $(ROOT_DIR)/components/lfric-xios

export UMPHYSICS_TARGETS = extract-um-physics
export UMPHYSICS_CLEAN   = clean-um-physics
export UM_FCM=$(UM_PHYS_DIR)/
export LFRIC_TARGET_PLATFORM ?= generic-x86
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN       ?= HEAD
META_FILE_DIR  = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build test-suite
	$(Q)echo > /dev/null

.PHONY: documentation doc docs
documentation doc docs: document-uml document-latex document-api
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
include build/project.mk

##############################################################################
# Launch test suite on XC40 ONLY
# temporary solution untill we can run this configuration on spice
#
# SUITE_CONFIG    - Path to rose-stem directory.
# SUITE_BASE_NAME - Name for suites.
# SUITE_GROUP_NAME_ABRV - Name(s) of the rose stem group(s) with abbreviations applied.
#
.PHONY: launch-test-suite-xc40
SUITE_GROUP ?= developer
ifneq (,$(findstring $(COMMA),$(SUITE_GROUP)))
  # Default for multiple groups: abbreviate names of groups in suite name
  SUITE_GROUP_ABRV ?= 1
else
  # Default for single group: keep full name of group in suite name
  SUITE_GROUP_ABRV ?= 0
endif
SUITE_GROUP_NAME_ABRV := $(subst $(COMMA),$(SPACE),$(shell echo $(SUITE_GROUP) | tr '[:lower:]' '[:upper:]'))
SUITE_GROUP_NAME_ABRV := $(subst $(SPACE),+,$(sort $(SUITE_GROUP_NAME_ABRV)))
SUITE_GROUP_NAME_ABRV := $(subst WEEKLY,W,$(SUITE_GROUP_NAME_ABRV))
SUITE_GROUP_NAME_ABRV := $(subst NIGHTLY-FULL-DEBUG,NFD,$(SUITE_GROUP_NAME_ABRV))
SUITE_GROUP_NAME_ABRV := $(subst NIGHTLY,N,$(SUITE_GROUP_NAME_ABRV))
SUITE_GROUP_NAME_ABRV := $(subst DEVELOPER,D,$(SUITE_GROUP_NAME_ABRV))
ifeq ($(SUITE_GROUP_ABRV),0)
launch-test-suite-xc40: SUITE_NAME = $(SUITE_BASE_NAME)-$$target-$(SUITE_GROUP)
else
launch-test-suite-xc40: SUITE_NAME = $(SUITE_BASE_NAME)-$$target-$(SUITE_GROUP_NAME_ABRV)
endif
ifdef VERBOSE
launch-test-suite-xc40: VERBOSE_ARG = --define-suite=VERBOSE=$(VERBOSE)
endif

launch-test-suite-xc40:
	$(Q)umask 022; for target in $(TEST_SUITE_TARGETS) ; do \
	if [[ $$target == "meto-xc40" || $$target == "meto-xcs" ]]; then\
           echo Launching $(PROJECT_NAME) test suite against $$target with $(SUITE_GROUP) group ; \
           rose stem --name=$(SUITE_NAME) \
                  --config=$(SUITE_CONFIG) \
                  --opt-conf-key=$$target --no-gcontrol \
                  $(CLEAN_OPT) $(QUIET_ARG) \
                  --define-suite=RDEF_PRECISION=$(RDEF_PRECISION) \
                  $(VERBOSE_ARG) \
                  --group=$(SUITE_GROUP); \
	else\
           echo Coupled configuration does NOT work for $$target target; \
        fi \
	done

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
test-suite: launch-suite-gscan launch-test-suite-xc40
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
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export CXX_LINK     = TRUE
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -print)))
build: export SCRATCH_DIR := $(WORKING_DIR)/um-scratch
build: export PROJECT      = $(PROJECT_NAME)
build: export WORKING_DIR := $(WORKING_DIR)/$(PROJECT_NAME)
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "fast-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS SAFE_OPTIMISATION FORTRAN_STANDARD
else ifeq "$(PROFILE)" "production"
build: export FFLAG_GROUPS = DEBUG WARNINGS RISKY_OPTIMISATION
else
  $(error Unrecognised profile "$(PROFILE)". Must be one of full-debug, fast-debug or production)
endif
build: ALWAYS
	$(call MESSAGE,========================================)
	$(call MESSAGE,Importing internal dependencies...)
	$(call MESSAGE,========================================)
	$(Q)for SUBPROJECT in $(INTERNAL_DEPENDENCIES) ; do \
		$(MAKE) $(QUIET_ARG) -f $$SUBPROJECT/build/import.mk ; done
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting coupler-oasis component)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
		SOURCE_DIR=$(ROOT_DIR)/components/coupler-oasis/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	         SOURCE_DIR=$(ROOT_DIR)/gungho/source \
	         WITHOUT_PROGRAMS=1
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting UM physics )
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(PROJECT_DIR)/build/extract_physics.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,Extracting $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(PROJECT_DIR)/build/extract_physics.mk
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk \
                  SOURCE_DIR=source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning coupler-oasis component)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	        SOURCE_DIR=$(ROOT_DIR)/components/coupler-oasis/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning Gungho dynamical core)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	        SOURCE_DIR=$(ROOT_DIR)/gungho/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning UM Physics)
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/um_physics/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/socrates/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(ROOT_DIR)/jules/source
	$(call MESSAGE,========================================)
	$(call MESSAGE,PSycloning $(PROJECT_NAME))
	$(call MESSAGE,========================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
                  SOURCE_DIR=source
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Generating $(PROJECT) namelist loaders)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/configuration.mk \
                  SOURCE_DIR=$(ROOT_DIR)/gungho/source \
                  META_FILE_DIR=$(META_FILE_DIR)
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Analysing $(PROJECT) build dependencies)
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(call MESSAGE,Extracting,"compile.mk", $(LFRIC_BUILD))
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk
	$(call MESSAGE,Extracting,"after compile.mk")

##############################################################################
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION                   \
                  = $(ROOT_DIR)/infrastructure/source        \
                    $(ROOT_DIR)/components/driver/source     \
                    $(ROOT_DIR)/components/science/source    \
                    $(ROOT_DIR)/components/inventory/source  \
                    $(ROOT_DIR)/gungho/source                \
                    $(ROOT_DIR)/um_physics/source            \
                    $(ROOT_DIR)/gungho/unit-test/support     \
                    $(ROOT_DIR)/components/lfric-xios/source \
                    $(ROOT_DIR)/components/coupler-oasis/source
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export META_FILE_DIR = rose-meta/lfric-lfric_coupled/HEAD
unit-tests/%: export PROGRAMS = $(PROJECT_NAME)_unit_tests
unit-tests/%: export PROJECT = $(PROJECT_NAME)
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-test
unit-tests: unit-tests/run


##############################################################################
# Integration tests
#
integration-tests:
	$(call MESSAGE,Testing,"There are no integration tests.")


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
