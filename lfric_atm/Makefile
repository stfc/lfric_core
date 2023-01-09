##############################################################################
# (c) Crown copyright 2017 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

PROJECT_NAME = lfric_atm

export PROFILE ?= fast-debug
export LFRIC_TARGET_PLATFORM ?= generic-x86

# This top level makefile is very order sensitive. Source code extraction and
# generation must happen in a certain order. Due to this we turn off
# multithreading for this file only. Any called recursively (i.e. with $(MAKE))
# run in parallel. Unless they specify NOTPARALLEL as well.
#
.NOTPARALLEL:

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..

export INTERNAL_DEPENDENCIES = $(ROOT_DIR)/infrastructure \
                               $(ROOT_DIR)/components/driver \
                               $(ROOT_DIR)/components/science \
                               $(ROOT_DIR)/components/inventory \
                               $(ROOT_DIR)/components/lfric-xios

export UMPHYSICS_TARGETS = extract-um-physics
export UMPHYSICS_CLEAN   = clean-um-physics
export UM_FCM=$(UM_PHYS_DIR)/
export OPTIMISATION_PATH ?= $(wildcard optimisation/$(LFRIC_TARGET_PLATFORM))

export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

META_VN             ?= HEAD
export META_FILE_DIR = $(PROJECT_DIR)/rose-meta/lfric-$(PROJECT_NAME)/$(META_VN)

.PHONY: default
default: build
	$(Q)echo > /dev/null

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)
include build/project.mk

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
# Test
#
.PHONY: test-suite
test-suite: SUITE_BASE_NAME = $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)
test-suite: SUITE_CONFIG = rose-stem
test-suite: launch-suite-gscan launch-test-suite
	$(Q)echo > /dev/null

##############################################################################
# Build
#
.PHONY: build
build: export BIN_DIR     ?= $(PROJECT_DIR)/bin
build: export COMPILE_OPTIONS = $(PWD)/build/compiler_options.mk
build: export CXX_LINK     = YES
build: export PROGRAMS    := $(basename $(notdir $(shell find source -maxdepth 1 -name '*.[Ff]90' -exec egrep -l "^\s*program" {} \;)))
build: export PROJECT      = lfric_atm
build: export SCRATCH_DIR := $(WORKING_DIR)/um-scratch
build: export WORKING_DIR := $(WORKING_DIR)/lfric_atm
ifeq "$(PROFILE)" "full-debug"
build: export FFLAG_GROUPS = DEBUG WARNINGS INIT RUNTIME NO_OPTIMISATION
else ifeq "$(PROFILE)" "fast-debug"
# Use fast-debug environment vars from build/fortran/*.mk configs
build: export FFLAG_GROUPS = DEBUG WARNINGS FASTD_INIT FASTD_RUNTIME SAFE_OPTIMISATION
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
	$(call MESSAGE,=========================================================)
	$(call MESSAGE,Compiling $(PROJECT))
	$(call MESSAGE,=========================================================)
	$Q$(MAKE) $(QUIET_ARG) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk

##############################################################################
# Unit tests
#
unit-tests:
	$(call MESSAGE,Testing,"There are no unit tests.")


##############################################################################
# Integration tests
#
integration-tests:
	$(call MESSAGE,Testing,"There are no integration tests.")


##############################################################################
# Clean
#
.PHONY: clean
clean: partial-clean
	$(call MESSAGE,Removing,"UM work space")
	$(Q)-rm -r $(WORKING_DIR)

.PHONY: partial-clean pclean
partial-clean pclean: ALWAYS
	$(call MESSAGE,Removing,"LFRic atmosphere  workspace")
	$(Q)-rm -r $(WORKING_DIR)/lfric_atm
	$(call MESSAGE,Removing,"LFRic atmosphere binaries")
	$(Q)-if [ -d bin ] ; then rm -r bin ; fi

