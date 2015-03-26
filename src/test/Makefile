##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

IGNORE_DEPENDENCIES = netcdf ESMF
EXTERNAL_LIBRARIES = netcdff esmf

export ROOT = ../..

include $(ROOT)/include.mk

PFUNIT_INSTALL_DIR = $(BUILD_DIR)/pfunit-install

SPECIAL_FFLAGS := $(filter-out -Werror,$(FFLAGS)) # pFUnit isn't that well put together

OBJ_DIR = $(BUILD_DIR)/tests
OBJ_SUBDIRS = $(patsubst %,$(OBJ_DIR)/%,$(filter-out data,$(shell find * -type d)))

DRIVER_OBJ=$(OBJ_DIR)/driver.o

TEST_SOURCES = $(shell find -name "*.pf")
TEST_OBJECTS = $(patsubst %.pf,$(OBJ_DIR)/%.o,$(TEST_SOURCES))

SUPPORT_SOURCE  = $(shell find -name "*.[Ff]90")
SUPPORT_OBJECTS = $(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(basename $(SUPPORT_SOURCE))))

SUBJECT_DIR     = $(BUILD_DIR)/dynamo
SUBJECT_SUBDIRS = $(shell find $(SUBJECT_DIR) -type d)
SUBJECT_MODINC  = $(patsubst %,-I %,$(SUBJECT_SUBDIRS))

.PHONY: test
test: $(OBJ_DIR)/test
	$(Q)cd $(OBJ_DIR); ./test -robust

.PRECIOUS: $(OBJ_DIR)/%.F90
$(OBJ_DIR)/%.F90: %.pf | $(OBJ_DIR) $(OBJ_SUBDIRS) $(DRIVER_OBJ)
	@echo Generating $@
	$(Q)$(PFUNIT_INSTALL_DIR)/bin/pFUnitParser.py $< $@

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.F90
	@echo Compiling $<
	$(Q)$(FC) $(FPPFLAGS) $(SPECIAL_FFLAGS) $(FFLAGS_TESTS) \
	          $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	          $(SUBJECT_MODINC) -I$(PFUNIT_INSTALL_DIR)/mod \
	          -c -o $@ $<

$(OBJ_DIR)/test: $(DRIVER_OBJ) $(SUPPORT_OBJECTS) $(TEST_OBJECTS)
	@echo Linking $@
	$(Q)$(FLINK) $(FLDFLAGS) -o $@ \
	             $^ \
	             $(SUBJECT_DIR)/modules.a \
	             $(PFUNIT_INSTALL_DIR)/lib/libpfunit.a \
	             $(patsubst %,-l%,$(EXTERNAL_LIBRARIES))

# Replace OBJ_SUBDIR with mangling of % ?
$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)/data $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	      $(SUBJECT_MODINC) -c -o $@ $<

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_SUBDIRS)
	$(Q)$(FC) $(FPPFLAGS) $(F_MOD_DESTINATION_ARG)$(OBJ_DIR)/$(dir $@) \
	      $(SUBJECT_MODINC) -c -o $@ $<

$(DRIVER_OBJ): | $(dir $(DRIVER_OBJ))
	$(MAKE) -f pfunit.mk PFUNIT_INSTALL_DIR=$(PFUNIT_INSTALL_DIR) \
	                     DRIVER_OBJ=$@

$(OBJ_DIR)/data: data
	@echo Copying data directory
	$(Q)cp -r $< $@

$(OBJ_DIR) $(OBJ_SUBDIRS):
	@echo "Creating $@"
	$(Q)mkdir -p $@

.PHONY: ALWAYS
ALWAYS:

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)
ifdef ALL
	$(MAKE) -f pfunit.mk clean
endif
