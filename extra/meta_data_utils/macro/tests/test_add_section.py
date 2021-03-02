##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
This module tests the macros found in add_section.py

Note - Rose is written using Python 2. These tests can be run from a terminal
by navigating into the macro directory and using:

python2 -m unittest tests.test_add_section
"""
import os
import sys
import unittest

# Add Rose library to path
sys.path.append("/home/h03/fcm/rose/lib/python/")
from add_section import AddField, AddStream
from rose.config import ConfigNode, load

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_CONF = "/rose-app.conf"
ADD_FIELD_EMPTY_CONF = "/add-field-empty-rose-app.conf"
ADD_FIELD_POPULATED_CONF = "/add-field-populated-rose-app.conf"
ADD_STREAM_POPULATED_CONF = "/add-stream-populated-rose-app.conf"
ADD_STREAM_EMPTY_CONF = "/add-stream-empty-rose-app.conf"


class TestAddSection(unittest.TestCase):
    def test_add_field_empty_stream(self):
        """Test a field is correctly added to an empty output stream"""
        old_conf = load(TEST_DIR + INPUT_CONF)
        expected_conf = load(TEST_DIR + ADD_FIELD_EMPTY_CONF)
        add_field = AddField().transform

        # Add field to output_stream(1)
        new_conf, reports = add_field(old_conf, output_stream=1)

        # Check the config is as expected
        self.assertEqual(expected_conf, new_conf)

        # Check number of and contents of reports
        self.assertEqual(len(reports), 1)
        self.assertEqual(reports[0].info, "Added field(1) to output_stream(1)")

    def test_add_field_no_stream(self):
        """Test that you cannot add a field to a nonexistent output stream"""
        old_conf = load(TEST_DIR + INPUT_CONF)
        add_field = AddField().transform

        with self.assertRaises(OSError):
            add_field(old_conf, output_stream=1.0)
        with self.assertRaises(OSError):
            add_field(old_conf, output_stream="one")
        with self.assertRaises(OSError):
            add_field(old_conf, output_stream=42)

    def test_add_field_populated_stream(self):
        """Test fields are correctly added to a stream with fields already
        present"""
        old_conf = load(TEST_DIR + INPUT_CONF)
        expected_conf = load(TEST_DIR + ADD_FIELD_POPULATED_CONF)
        add_field = AddField().transform

        new_conf, reports = add_field(old_conf, output_stream=2)

        self.assertEqual(expected_conf, new_conf)

        self.assertEqual(len(reports), 1)
        self.assertEqual(reports[0].info, "Added field(2) to output_stream(2)")

    def test_add_stream_empty(self):
        """Test that a new stream is correctly added to a configuration that
        does not already contain any streams"""
        old_conf = ConfigNode()
        expected_conf = load(TEST_DIR + ADD_STREAM_EMPTY_CONF)
        add_stream = AddStream().transform

        new_conf, reports = add_stream(old_conf)

        self.assertEqual(expected_conf, new_conf)

        self.assertEqual(len(reports), 2)
        self.assertEqual(reports[0].info, "Added output_stream(1)")
        self.assertEqual(reports[1].info, "Added field(1) to output_stream(1)")

    def test_add_stream_populated(self):
        """Test stream is added correctly when there are already stream
         present"""
        old_conf = load(TEST_DIR + INPUT_CONF)
        expected_conf = load(TEST_DIR + ADD_STREAM_POPULATED_CONF)
        add_stream = AddStream().transform

        new_conf, reports = add_stream(old_conf)

        self.assertEqual(expected_conf, new_conf)

        self.assertEqual(len(reports), 2)
        self.assertEqual(reports[0].info, "Added output_stream(3)")
        self.assertEqual(reports[1].info, "Added field(1) to output_stream(3)")
