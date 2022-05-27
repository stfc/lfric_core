##############################################################################
# (c) Crown copyright 2021 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic um_physics project.
# Targets provided are detailed below.
# TODO: Currently all targets are interactive only, #2757 will add the testing
#		to a rose-stem test-suite.
#
# all: (default) Complete build and test the um_physics.
#      Intended for interactive use. Currently only builds and runs unit-tests.
# unit-tests: Run unit tests interactively.
# clean: Delete all final products and working files.
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
#
##############################################################################

PROJECT_NAME = um_physics

PROFILE ?= fast-debug

export PROJECT_DIR := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))
export ROOT_DIR    ?= ..

export INTERNAL_DEPENDENCIES = $(ROOT_DIR)/infrastructure

.PHONY: all
all: default

.PHONY: default
default: unit-tests
	$(Q)echo > /dev/null

META_VN ?= HEAD

include $(ROOT_DIR)/infrastructure/build/lfric.mk
include $(INTERNAL_DEPENDENCIES:=/build/import.mk)


##############################################################################
# Unit tests
#
unit-tests/%: export ADDITIONAL_EXTRACTION                   \
                  = $(ROOT_DIR)/infrastructure/source        \
                    $(ROOT_DIR)/components/driver/source     \
                    $(ROOT_DIR)/components/lfric-xios/source \
                    $(ROOT_DIR)/components/science/source
unit-tests/%: export BIN_DIR ?= $(PROJECT_DIR)/test
unit-tests/%: export EXTERNAL_STATIC_LIBRARIES += pfunit
unit-tests/%: export PROGRAMS = um_physics_unit_tests
unit-tests/%: export PROJECT = um_physics
unit-tests/%: export SOURCE_DIR = source
unit-tests/%: export TEST_DIR = unit-test
unit-tests/%: export WORKING_DIR := $(WORKING_DIR)/unit-tests
unit-tests: unit-tests/run

##############################################################################
# Create a consolidated rose-meta file for um_physics
#
create_rose-meta: rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf
rose-meta/$(PROJECT_NAME)/HEAD/rose-meta.conf:
	$(call MESSAGE, Creating, $@)
	$(Q)mkdir -p $(dir $@)
	$(Q)for meta_file in $(addprefix $(PROJECT_NAME)/,$(addsuffix $(META_VN), \
	$(dir $(wildcard rose-meta/*/.)))); \
	do echo import\=$${meta_file} >> $@; done

##############################################################################
# Clean
#
.PHONY: clean
clean: ALWAYS
	$(call MESSAGE,Removing,"um_physics work space")
	$(Q)-rm -r $(WORKING_DIR)
	$(call MESSAGE,Removing,"um_physics generated rose-meta file")
	$(Q)-rm -r rose-meta/$(PROJECT_NAME)
	$(call MESSAGE,Removing,"um_physics binaries")
	$(Q)-if [ -d test ] ; then rm -r test ; fi
