##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains macros to add output streams and their diagnostic
fields to LFRic configurations in Rose config editor"""
import rose.macro


class AddField(rose.macro.MacroBase):
    """Adds a field to the given output_stream"""

    def transform(self, config, meta_config=None, output_stream=1):
        """Adds a field to the given output_stream"""

        if config.get(["output_stream({})".format(output_stream)]) is None:
            raise OSError("output_stream({}) not found".format(output_stream))
        field = 1
        added = False

        while not added:
            section_name = "output_stream({}):field({})".format(output_stream,
                                                                field)
            if config.get([section_name]) is None:
                config.set([section_name])
                self.add_report(section_name, None, None,
                                "Added field({}) to output_stream({})".format(
                                    field, output_stream))
                added = True
            field += 1
        return config, self.reports


class AddStream(rose.macro.MacroBase):
    """Adds a new output_stream"""

    def transform(self, config, meta_config=None):
        """Adds a new output_stream"""

        output_stream = 1
        added = False

        while not added:
            section_name = "output_stream({})".format(output_stream)
            subsection_name = "output_stream({}):field(1)".\
                format(output_stream)

            if config.get([section_name]) is None:
                config.set([section_name])
                self.add_report(section_name, None, None,
                                "Added output_stream({})".
                                format(output_stream))
                config.set([subsection_name])
                self.add_report(subsection_name, None, None,
                                "Added field(1) to output_stream({})".format(
                                    output_stream))
                added = True
            output_stream += 1
        return config, self.reports
