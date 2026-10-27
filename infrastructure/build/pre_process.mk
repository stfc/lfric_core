##############################################################################
# Copyright (c) 2026,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE which you
# should have received as part of this distribution.
##############################################################################
#
# Locate and compile all Fortran .F90 and .X90 files and preprocess them.
#
# The following variables may be specified to modify the compile process...
#
# PRE_PROCESS_INCLUDE_DIRS: Space separated list of directories to search for
#                           inclusions.
# PRE_PROCESS_MACROS: Space separated list of macro definitions in the form
#                     NAME[=MACRO] to be passed to the compiler.


# ------
# ROOT: Project directory
# BIN_DIR: Path to directory for resulting executables.
#          Default: $(ROOT)/bin
# FFLAG_GROUPS: Space separated list of FFLAG_<group name> variables to use in
#               building up the FFLAGS variable. Only the group name is
#               specified.
# LINK_TYPE: 'static' or 'dynamic' linking.
#            Default: dynamic, except on Crays where it's static
# PROGRAMS: Names of programs to compile.
#           Default: Everything listed in programs.mk
# PROJECT_MAKE_DIR: Used to locate project specific targets modifiers and
#                   such.
# COMPILE_OPTIONS: Name of an optional file that can be included to list
#                  project specific compile options
#
##############################################################################

.SECONDEXPANSION:

# Build a set of "-I" arguments to seach the whole object tree:
INCLUDE_ARGS := $(subst ./,-I,$(shell find . -mindepth 1 -type d -print)) \
                $(addprefix -I, $(PRE_PROCESS_INCLUDE_DIRS))

# Build a set of "-D" argument for any pre-processor macros
#
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

PP_FILES := $(shell find $(WORKING_DIR) -name '*.[FX]90' -print)
# Convert to lower case
LOWER_PP_FILES := $(patsubst %.X90,%.x90,$(patsubst %.F90,%.f90,$(PP_FILES)))

pre_process: $(LOWER_PP_FILES)

FPPFLAGS += -P

# nvfortran as preprocessor needs different options
PREPROCESSOR = $(firstword $(subst -, ,$(notdir $(FPP))))

# Debug output
$(info FC $(FC) PREPROCESSOR $(PREPROCESSOR) FPP $(FPP) FPPFLAGS $(FPPFLAGS))

ifeq ("$(PREPROCESSOR)", "nvfortran")
%.x90: /%.X90
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) -o $@ $<
	$Qrm $<
%.f90: %.F90
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) -o $@ $<
	$Qrm $<
else
%.x90: %.X90
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) $< $@
	$Qrm $<
%.f90: %.F90
	$(call MESSAGE,Preprocessing, $<)
	$Q$(FPP) $(FPPFLAGS) $(MACRO_ARGS) $< $@
	$Qrm $<
endif

# Get MESSAGE function from lfric.mk
include $(LFRIC_BUILD)/lfric.mk
