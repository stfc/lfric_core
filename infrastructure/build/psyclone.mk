##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
#
# Run this make file to generate PsyKAl source in WORKING_DIR from algorithms
# and kernels in SOURCE_DIR.
#
.PHONY: psykal-source
psykal-source: jump

include $(LFRIC_BUILD)/lfric.mk
include $(LFRIC_BUILD)/fortran.mk

ALGORITHM_ONLY_FILES := ${shell find $(SOURCE_DIR) -name '*.x90' -exec grep -EiL '[^!]\s*call\s+invoke\s*[$(OPEN_PAREN)]' {} \; }
HAND_WRITTEN_FILES := $(shell find $(SOURCE_DIR) -path '*/psy/*' -a -name '*_psy.[Ff]90' -print)
OVERRIDDEN_FILES = $(patsubst %_psy.f90,%.x90,$(subst /psy/,/algorithm/,$(HAND_WRITTEN_FILES)))
ALGORITHM_ONLY_SOURCE = $(patsubst $(SOURCE_DIR)/%.x90,$(WORKING_DIR)/%.f90,$(ALGORITHM_ONLY_FILES) $(OVERRIDDEN_FILES))
PSY_SOURCE := $(patsubst $(SOURCE_DIR)/%.x90,$(WORKING_DIR)/%_psy.f90,$(filter-out $(ALGORITHM_ONLY_FILES) $(OVERRIDDEN_FILES),$(shell find $(SOURCE_DIR) -name '*.x90')))
KERNEL_SOURCE := $(patsubst ./%,$(WORKING_DIR)/%,$(shell find $(SOURCE_DIR) -path '*/kernel/*' -a -name '*.[Ff]90' -print))
MACRO_ARGS := $(addprefix -D,$(PRE_PROCESS_MACROS))

DIRECTORIES := $(patsubst $(SOURCE_DIR)/%,$(WORKING_DIR)/%,$(shell find $(SOURCE_DIR) -type d -print))

.PHONY: jump
jump: $(PSY_SOURCE) $(ALGORITHM_ONLY_SOURCE)

ifdef OPTIMISATION_PATH
  GLOBAL_OPTIMISATION_FILE = $(OPTIMISATION_PATH)/global.py
  LOCAL_OPTIMISATION_FILE = $(OPTIMISATION_PATH)/$*.py

$(WORKING_DIR)/%_psy.f90: $(WORKING_DIR)/%.x90 $$(LOCAL_OPTIMISATION_FILE) \
                          | $(KERNEL_SOURCE)
	$(call MESSAGE,Full psyclone - local optimisations,$<)
	$(Q)psyclone -api dynamo0.3 -l all -d $(WORKING_DIR) \
	             -s $(LOCAL_OPTIMISATION_FILE) \
	             -opsy $@ -oalg $(WORKING_DIR)/$*.f90 $<

$(WORKING_DIR)/%_psy.f90: $(WORKING_DIR)/%.x90 $(GLOBAL_OPTIMISATION_FILE) \
                          | $(KERNEL_SOURCE)
	$(call MESSAGE,Full psyclone - global optimisations,$<)
	$(Q)psyclone -api dynamo0.3 -l all -d $(WORKING_DIR) \
	             -s $(GLOBAL_OPTIMISATION_FILE) \
	             -opsy $@ -oalg $(WORKING_DIR)/$*.f90 $<

$(OPTIMISATION_PATH)/global.py:
	$(error File not found: $(OPTIMISATION_PATH)/global.py)

else
$(WORKING_DIR)/%_psy.f90: $(WORKING_DIR)/%.x90 | $(KERNEL_SOURCE)
	$(call MESSAGE,Full psyclone,$(subst $(PWD)/,,$<))
	$(Q)psyclone -api dynamo0.3 -l all -d $(WORKING_DIR) \
	             -opsy $@ -oalg $(WORKING_DIR)/$*.f90 $<
endif

$(WORKING_DIR)/%.f90: $(WORKING_DIR)/%.x90
	$(call MESSAGE,Algorithm only,$<)
	$(Q)psyclone -api dynamo0.3 -l all -d $(WORKING_DIR) \
	             -oalg $@ -opsy /dev/null $<

.PRECIOUS: $(WORKING_DIR)/%.x90
$(WORKING_DIR)/%.x90: $(SOURCE_DIR)/%.x90 \
                      | $$(patsubst $$(PERCENT)/,$$(PERCENT),$$(dir $$@))
	$(call MESSAGE,Preprocessing,$<)
	$(Q)$(FPP) $(FPPFLAGS) $(MACRO_ARGS) $< $@

$(DIRECTORIES): ALWAYS
	$(call MESSAGE,Creating,$@)
	$(Q)mkdir -p $@
