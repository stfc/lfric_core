##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
export IGNORE_DEPENDENCIES += netcdf MPI yaxt pfunit_mod mod_oasis
export EXTERNAL_DYNAMIC_LIBRARIES += yaxt yaxt_c netcdff netcdf hdf5

TEMPLATE_TOOL = $(LFRIC_BUILD)/tools/Templaterator
TYPE_TABLE_real32 = real
TYPE_TABLE_real64 = real
TYPE_TABLE_int32 = integer

.PHONY: import-infrastructure
import-infrastructure: $(WORKING_DIR)/field/field_real32_mod.f90 \
                       $(WORKING_DIR)/field/field_real64_mod.f90 \
		       $(WORKING_DIR)/field/field_int32_mod.f90
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/extract.mk SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source
	$Q$(MAKE) $(QUIET_ARG) -f $(LFRIC_BUILD)/psyclone/psyclone.mk \
	          SOURCE_DIR=$(LFRIC_INFRASTRUCTURE)/source \
	          OPTIMISATION_PATH=$(OPTIMISATION_PATH)

$(WORKING_DIR)/field/field_%_mod.f90: $(LFRIC_INFRASTRUCTURE)/source/field/field_mod.t90 \
                                      | $(WORKING_DIR)/field
	$(call MESSAGE, Templating, $<)
	$Q$(TEMPLATE_TOOL) $< -o $@ -s type=$(TYPE_TABLE_$*) -s kind=$*

$(WORKING_DIR)/field:
	$Qmkdir -p $@
