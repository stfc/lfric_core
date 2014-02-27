ROOT = ../..
include $(ROOT)/Makefile.inc

OBJ_DIR    = $(BUILD_DIR)/dynamo
TEST_SRC   = $(wildcard *.pf)
TEST_OBJS  = $(patsubst %.pf,$(OBJ_DIR)/%.o,$(TEST_SRC))
SRC        = $(wildcard *.[Ff]90)
OBJS       = $(patsubst %.[Ff]90,$(OBJ_DIR)/%.o,$(SRC))
DRIVER_OBJ = $(OBJ_DIR)/driver.o
TARGET_OBJ = $(patsubst %Test.pf,$(OBJ_DIR)/%.o,$(TEST_SRC))

.PHONY: test
test: pfunit $(OBJ_DIR)/test
	$(OBJ_DIR)/test

include pfunit.mk

$(OBJ_DIR)/test: $(TEST_OBJS) $(OBJS) $(DRIVER_OBJ) $(TARGET_OBJ)
	$(FC03) -L$(PFUNIT_INSTALL)/lib -o $@ $^ -lpfunit

$(DRIVER_OBJ): $(PFUNIT_INSTALL)/include/driver.F90 testSuites.inc
	$(FC03) $(F03FLAGS) -c -I$(PFUNIT_INSTALL)/include -I. \
	        $(F03MODULES) $(PFUNIT_INSTALL)/mod -o $@ $<

$(OBJ_DIR)/%.o: $(OBJ_DIR)/%.F90
	$(FC03) $(F03FLAGS) -c -I$(PFUNIT_INSTALL)/mod \
	        $(F03MODULES) $(OBJ_DIR) -o $@ $<

.PRECIOUS: $(OBJ_DIR)/%.F90
$(OBJ_DIR)/%.F90: %.pf
	$(PFUNIT_INSTALL)/bin/pFUnitParser.py $< $@

.PHONY: clean
clean: cleanpfunit
	-rm $(OBJS)
	-rm $(TEST_OBJS)
	-rm $(OBJ_DIR)/test
