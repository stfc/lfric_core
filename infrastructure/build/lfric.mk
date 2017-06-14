##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Include this file from your model make file in order to gain access to the
# LFRic build system. Include it at the end of the make file as it contains
# targets which you do not want to become the default target.
#
# Macros provided by including this file...
#
# LFRIC_BUILD: Path to the build system
#
# Macros expected by the build system are as follows...
#
# ROOT: Absolute path to the project directory
# OPTIMISATION_PATH: Where PSyclone optimisation scripts may be found.
# WORKING_DIR: Path to scratch space in which intermediate files will be
#              placed. This should be somewhere with good "many small
#              files" performance, i.e. probably not Lustre.
#              Default: ./working
# VERBOSE: Set in order to see actual commands issued by the build system.
#
# Plus the normal compiler macros...
#
# FC: Fortran compiler command
# FFLAGS: Fortran compiler flags including "-I" module file locations
# LD: Linker command, probably the same as FC
# MPILD: MPI linker command
# LDFLAGS: Linker flags including "-L" library locations
#
##############################################################################

.SECONDEXPANSION:

# Ensure make offers the features we need...
#
$(info ** Make version $(MAKE_VERSION))
ifeq ($(filter else-if,$(value .FEATURES)),)
  $(error The build system requires else-if support from GMake)
endif

# The default value of FC is almost always "f77" which is of no use to us.
# An empty FC is also of no use.
ifneq "$(or $(filter default, $(origin FC)), $(filter x, x$(FC)))" ""
  $(error The FC environment variable must be set to a Fortran compiler command)
endif

# Default variables...
#
export WORKING_DIR ?= working

# Make the build system available...
#
export LFRIC_BUILD := $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

# Make compiler macros available...
#
# Sometimes FC holds a full path which needs to be stripped off. It may also
# include a version number which also needs to go.
#
FORTRAN_COMPILER := $(firstword $(subst -, ,$(notdir $(FC))))

# Attempt to identify Cray systems...
#
ifdef PE_ENV
  CRAY_ENVIRONMENT = true
  ifeq '$(PE_ENV)' 'CRAY'
    FORTRAN_COMPILER = crayftn
  else ifeq '$(PE_ENV)' 'INTEL'
    FORTRAN_COMPILER = ifort
  else ifeq '$(PE_ENV)' 'GNU'
    FORTRAN_COMPILER = gfortran
  else ifeq '$(PE_ENV)' 'PGI'
    FORTRAN_COMPILER = pgfortran
  else
    $(error Unrecognised Cray programming environment)
  endif
endif

include $(LFRIC_BUILD)/fortran/$(FORTRAN_COMPILER).mk
export F_MOD_DESTINATION_ARG OPENMP_ARG

FFLAGS += $(FFLAGS_COMPILER)
export FFLAGS

# Set up verbose logging...
#
ifdef VERBOSE
  Q :=
  VERBOSE_ARG = -verbose
  SHORT_VERBOSE_ARG = -v
  DOUBLE_VERBOSE_ARG = --verbose
else
  Q := @
  VERBOSE_REDIRECT = >/dev/null
endif
export Q

# We only want to send terminal control characters if there is a terminal to
# interpret them...
#
ifneq 'x$(TERM)' 'x'
  MESSAGE = @echo -e \\x1b[1m$(1)\\x1b[0m $(2)
else
  MESSAGE = @echo *$(1)* $(2)
endif

# Set up some special macros for hard to escape characters
#
EMPTY :=
SPACE := $(EMPTY) # This comment highlights space character.
PERCENT := %
OPEN_PAREN := (

# Prerequisite for targets which should always be run.
#
.PHONY: ALWAYS
ALWAYS:

# The directory containing the target. Useful for order-only prerequisites to
# create that directory.
#
TARGET_DIR = $(patsubst $(PERCENT)/,$(PERCENT),$(dir $@))

##############################################################################
# Build API documentation
#
PHONY: api-documentation-%
api-documentation-%:
	$(call MESSAGE,API,$*)
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)( cat $(CONFIG_DIR)/Doxyfile; \
	      echo INPUT = $(SOURCE_DIR); \
	      echo OUTPUT_DIRECTORY = $(DOCUMENT_DIR) ) \
	    | doxygen - $(VERBOSE_REDIRECT)

##############################################################################
# Build UML documentation
#
.PHONY: uml-documentation-%
uml-documentation-%: $$(patsubst $$(SOURCE_DIR)/$$(PERCENT).puml,$$(DOCUMENT_DIR)/$$(PERCENT).pdf,$$(wildcard $$(SOURCE_DIR)/*.puml))
	$(Q)echo >/dev/null

.PRECIOUS: $(DOCUMENT_DIR)/%.pdf
$(DOCUMENT_DIR)/%.pdf: $(DOCUMENT_DIR)/%.svg
	$(call MESSAGE,Translating,$@)
	$(Q)inkscape $< --export-pdf=$@

.PRECIOUS: $(DOCUMENT_DIR)/%.svg
$(DOCUMENT_DIR)/%.svg: $(SOURCE_DIR)/%.puml \
                      $$(addprefix $$(SOURCE_DIR)/,$$(shell sed -n -e 's/!include[ ]*\([^ \n]*\)/\1/p' $$(SOURCE_DIR)/$$*.puml))
	$(call MESSAGE,Generating,$@)
	$(Q)mkdir -p $(DOCUMENT_DIR)
	$(Q)plantuml $(SHORT_VERBOSE_ARG) -tsvg -o $(abspath $(dir $@)) $(abspath $<)

##############################################################################
# Run integration tests.
#
.PHONY: run-integration-test-%
run-integration-test-%: PYTHONPATH := $(PYTHONPATH):$(LFRIC_BUILD)
run-integration-test-%: $(patsubst %,run-run-integration-test-%,$(PROGRAMS))
	$(Q)echo >/dev/null

.PHONY: run-run-integration-test-%
run-run-integration-test-%: compile
	$(call MESSAGE,Running,$*)
	$(Q)$(patsubst %,$(SOURCE_DIR)/%.py,$*) $(patsubst %,$(BIN_DIR)/%,$*)

##############################################################################
# Run unit tests.
#
.PHONY: run-unit-test-%
run-unit-test-%: compile
	$(call MESSAGE,Running,$(PROGRAMS))
	$(Q)cd $(WORKING_DIR); mpiexec -n 1 $(BIN_DIR)/$(PROGRAMS) $(DOUBLE_VERBOSE_ARG)

##############################################################################
# Simple build process targets.
#

.PHONY: compile
compile:
	$(MAKE) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/analyse.mk
	$(MAKE) -C $(WORKING_DIR) -f $(LFRIC_BUILD)/compile.mk \
	        BIN_DIR=$(BIN_DIR) PROGRAMS="$(PROGRAMS)"

##############################################################################
# Generate configuration source.
#
.PHONY: configuration
configuration:
	$(MAKE) -f $(LFRIC_BUILD)/configuration.mk

##############################################################################
# Generate pFUnit unit tests.
#
.PHONY: pfunit
pfunit:
	$(MAKE) -f $(LFRIC_BUILD)/pfunit.mk

##############################################################################
# Generate PSyKAl source.
#
.PHONY: psykal
psykal:
	$(MAKE) -f $(LFRIC_BUILD)/psyclone.mk
