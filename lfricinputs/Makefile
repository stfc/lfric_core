##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
#
# Master make file for LFRic Inputs project.
# Targets provided our detailed below...
#
# test-suite: (default) Run rose stem using "developer" group of tasks

PROJECT_NAME = lfricinputs

export ROOT_DIR   ?= ..
export SUITE_GROUP ?= developer
export SUITE_GROUP_NAME ?= $(notdir $(realpath $(shell pwd)/$(ROOT_DIR)))-$(PROJECT_NAME)-.*

include $(ROOT_DIR)/infrastructure/build/lfric.mk

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
# Launch test suite
#
# SUITE_CONFIG    - Path to rose-stem directory.
# SUITE_BASE_NAME - Name for suites.
#
.PHONY: launch-test-suite
launch-test-suite: SUITE_GROUP ?= developer
launch-test-suite: SUITE_NAME = $(SUITE_BASE_NAME)-$(SUITE_GROUP)
ifdef VERBOSE
launch-test-suite: VERBOSE_ARG = --define-suite=VERBOSE=$(VERBOSE)
endif
launch-test-suite:
	$(Q)umask 022;  \
	echo Launching $(PROJECT_NAME) test suite with $(SUITE_GROUP) group ; \
	rose stem --name=$(SUITE_NAME) \
	          --config=$(SUITE_CONFIG) \
	          --no-gcontrol \
	          $(CLEAN_OPT) $(QUIET_ARG) \
	          --define-suite=RDEF_PRECISION=$(RDEF_PRECISION) \
                  $(VERBOSE_ARG) \
	          --group=$(SUITE_GROUP);

