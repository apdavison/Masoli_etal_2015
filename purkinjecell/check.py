"""
Check the ion channel distributions in the PyNN version of the model
by comparing them to those in the ModelDB version
"""


from Purkinje import Purkinje
from Purkinje2 import population

import os
import sys

from neuron import psection
from parameters import ParameterSet as P


class OutputGrabber(object):
    """
    Class used to grab standard output or another stream.

    Thanks to Devan Williams: https://stackoverflow.com/a/29834357
    """
    escape_char = "\b"

    def __init__(self, stream=None):
        self.origstream = stream
        if self.origstream is None:
            self.origstream = sys.stdout
        self.origstreamfd = self.origstream.fileno()
        self.capturedtext = ""
        # Create a pipe so the stream can be captured:
        self.pipe_out, self.pipe_in = os.pipe()

    def __enter__(self):
        self.start()
        return self

    def __exit__(self, type, value, traceback):
        self.stop()

    def start(self):
        """
        Start capturing the stream data.
        """
        self.capturedtext = ""
        # Save a copy of the stream:
        self.streamfd = os.dup(self.origstreamfd)
        # Replace the original stream with our write pipe:
        os.dup2(self.pipe_in, self.origstreamfd)

    def stop(self):
        """
        Stop capturing the stream data and save the text in `capturedtext`.
        """
        # Print the escape character to make the readOutput method stop:
        self.origstream.write(self.escape_char)
        # Flush the stream to make sure all our data goes in before
        # the escape character:
        self.origstream.flush()
        self.readOutput()
        # Close the pipe:
        os.close(self.pipe_out)
        # Restore the original stream:
        os.dup2(self.streamfd, self.origstreamfd)

    def readOutput(self):
        """
        Read the stream data (one byte at a time)
        and save the text in `capturedtext`.
        """
        while True:
            char = os.read(self.pipe_out, 1)
            if not char or self.escape_char in char:
                break
            self.capturedtext += char


def get_section_properties(sec):
    out = OutputGrabber(sys.stdout)
    out.start()
    psection(sec)
    out.stop()
    raw_output = out.capturedtext

    lines = raw_output.split("\n")
    name, params = lines[0].split(" { ")
    for s in params.split():
        exec(s)  # nseg, L, Ra
    output = {
        name: {
            "nseg": nseg,
            "L": L,
            "Ra": Ra,
            "mechanisms": {}
        }
    }
    i = 1
    for line in lines[1:]:
        if "insert" in line:
            # consider doing this with a regexp
            start = line.find("insert") + 7
            lbracket = line.find("{")
            rbracket = line.find("}")
            mech_name = line[start:lbracket - 1]
            output[name]["mechanisms"][mech_name] = {}
            for s in line[lbracket + 1:rbracket].split():
                pname, value = s.split("=")
                output[name]["mechanisms"][mech_name][pname] = float(value)
    return output


purkinje1 = Purkinje()
purkinje2 = population[0]._cell

soma_properties_orig = get_section_properties(purkinje1.soma)
soma_properties_new = get_section_properties(purkinje2.section_labels['soma'])

diff = P(soma_properties_new) - P(soma_properties_orig)
assert diff == ({}, {})

diffs = []

for dend1, dend2_index in zip(purkinje1.dend, purkinje2.morphology.section_groups['dend']):
    name = dend1.name()
    #print(name)
    dend2_id = purkinje2.morphology.segments[dend2_index].id  # !! indirect much?
    dend2 = purkinje2.sections[dend2_id]
    assert dend2.name() == name

    properties_orig = get_section_properties(dend1)
    properties_new = get_section_properties(dend2)
    # round up L to 3 decimal places to avoid spurious diffs from floating point errors
    properties_orig[name]['L'] = round(properties_orig[name]['L'], 3)
    properties_new[name]['L'] = round(properties_new[name]['L'], 3)

    diff = P(properties_new) - P(properties_orig)
    if diff != ({}, {}):
        diffs.append(diff)

# axon
for name in ("axonAIS", "axonAISK", "axonNOR", "axonNOR2", "axonNOR3", "axonmyelin",
             "axonmyelin2", "axonmyelin3", "axonmyelin4", "axoncoll", "axoncoll2"):
    axon1 = getattr(purkinje1, name)
    axon2 = purkinje2.section_labels[name]

    properties_orig = get_section_properties(axon1)
    properties_new = get_section_properties(axon2)

    diff = P(properties_new) - P(properties_orig)
    if diff != ({}, {}):
        diffs.append(diff)

assert len(diffs) == 0
