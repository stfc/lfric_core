##############################################################################
# (C) British Crown Copyright 2018 Met Office.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rose. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

export ROOT_DIR ?= ../../

.PHONY: test
test: export GPL_DIR ?= $(ROOT_DIR)/GPL
test: export GPL_PYTHON_LIB ?= $(GPL_DIR)/lib/python
test: PYTHONPATH:=$(PYTHONPATH):./:$(GPL_PYTHON_LIB)
test:
	python -m unittest discover -p '*test.py' --verbose
