##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

export ROOT_DIR   ?= ../../..
export PATH       := $(PATH):$(abspath $(ROOT_DIR)/GPL/bin)
export PYTHONPATH := $(PYTHONPATH):$(abspath $(ROOT_DIR)/GPL/lib/python)

.PHONY: test
test:
	python2.7 -m unittest discover -p '*test.py' --verbose
