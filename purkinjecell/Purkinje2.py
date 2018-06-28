"""
PyNN version of the Purkinje cell model from Masoli et al.,2015

Author: Andrew P. Davison, 2018
"""

import sys
import io
from collections import defaultdict
from neuron import h  # for debugging
import numpy as np
from neuroml import Morphology, Segment, Point3DWithDiam as P
from pyNN.morphology import NeuroMLMorphology, uniform, with_label, any
from pyNN.parameters import IonicSpecies
import pyNN.neuron as sim
from pyNN.neuron.nmodl import NMODLChannel

from PC_param import pc_param 
from ion_channel_params import ion_channel_parameters

use_arraymorph = False

soma = Segment(proximal=P(x=0, y=0, z=0, diameter=29.8),
               distal=P(x=0, y=29.8, z=0, diameter=29.8),
               name="soma")

dendrites = []
with open("PC_dendnames.dlist") as fp:
    dendrite_section_names = [line.strip() for line in fp if len(line) > 1]
section_coordinates = np.genfromtxt("coordinate.csv")
section_connections = np.genfromtxt("connections.csv")
assert len(dendrite_section_names) == section_coordinates.shape[0]

if use_arraymorph:
    raise NotImplementedError("to do")
else:
    for dend_name, C in zip(dendrite_section_names, section_coordinates):
        dendrites.append(
            Segment(proximal=P(x=C[1], y=C[2], z=C[3], diameter=C[4]),
                    distal=P(x=C[5], y=C[6], z=C[7], diameter=C[8]),
                    name=dend_name)
        )
    dendrites[0].parent = soma
    assert (section_connections[:, 1] == 0).all()
    assert (section_connections[:, 3] == 1).all()
    for c in section_connections:
        dendrites[int(c[0])].parent = dendrites[int(c[2])]


axonAIS = Segment(proximal=P(x=0, y=0, z=0, diameter=0.97),
                  distal=P(x=17, y=0, z=0, diameter=0.97),
                  name='axonAIS', parent=soma)
axonAISK = Segment(proximal=P(x=17, y=0, z=0, diameter=0.97),
                   distal=P(x=21, y=0, z=0, diameter=0.97),
                   name='axonAISK', parent=axonAIS)
axonmyelin1 = Segment(proximal=P(x=21, y=0, z=0, diameter=0.73),
                   distal=P(x=121, y=0, z=0, diameter=0.73),
                   name='axonmyelin', parent=axonAISK)
axonNOR1 = Segment(proximal=P(x=121, y=0, z=0, diameter=0.73),
                   distal=P(x=125, y=0, z=0, diameter=0.73),
                   name='axonNOR', parent=axonmyelin1)
axonmyelin2 = Segment(proximal=P(x=125, y=0, z=0, diameter=0.73),
                   distal=P(x=225, y=0, z=0, diameter=0.73),
                   name='axonmyelin2', parent=axonNOR1)
axonNOR2 = Segment(proximal=P(x=225, y=0, z=0, diameter=0.73),
                   distal=P(x=229, y=0, z=0, diameter=0.73),
                   name='axonNOR2', parent=axonmyelin2)
axonmyelin3 = Segment(proximal=P(x=229, y=0, z=0, diameter=0.73),
                   distal=P(x=329, y=0, z=0, diameter=0.73),
                   name='axonmyelin3', parent=axonNOR2)
axonNOR3 = Segment(proximal=P(x=329, y=0, z=0, diameter=0.73),
                   distal=P(x=333, y=0, z=0, diameter=0.73),
                   name='axonNOR3', parent=axonmyelin3)
axonmyelin4 = Segment(proximal=P(x=333, y=0, z=0, diameter=0.73),
                   distal=P(x=433, y=0, z=0, diameter=0.73),
                   name='axonmyelin4', parent=axonNOR3)
axoncoll1 = Segment(proximal=P(x=229, y=0, z=0, diameter=0.60),
                   distal=P(x=229, y=0, z=100, diameter=0.60),
                   name='axoncoll', parent=axonNOR2)
axoncoll2 = Segment(proximal=P(x=229, y=0, z=100, diameter=0.60),
                   distal=P(x=229, y=0, z=200, diameter=0.60),
                   name='axoncoll2', parent=axoncoll1)
axon = [axonAIS, axonAISK, axonmyelin1, axonNOR1, axonmyelin2, axonNOR2,
        axonmyelin3, axonNOR3, axonmyelin4, axoncoll1, axoncoll2]


segments = [soma] + dendrites + axon
for i, seg in enumerate(segments):
    seg.id = i
morph = NeuroMLMorphology(Morphology(segments=segments))
morph.section_groups["dend"] = np.arange(1, 1 + len(dendrites))  # dendrites

subsets = np.genfromtxt("ModelViewParmSubset.txt", dtype=int)
for subset_id in range(88):
    morph.section_groups["dend_subset{}".format(subset_id)] = \
        1 + subsets[:, 0][subsets[:, 1] == subset_id]  # offset is because dend indices start at 1

Purkinje = sim.MultiCompartmentNeuron  # standard base class for multi-compartment neurons
Purkinje.label = "PurkinjeNeuron"
Purkinje.ion_channels = {
    name: NMODLChannel(name)
    for name in ('Leak', 'HCN1', 'Nav1_6', 'Kv3_4', 'Kv1_1',
                'Cav3_2', 'Kca3_1', 'Kir2_3', 'Kca1_1', 'Kca2_2',
                'Kv4_3', 'Kv3_3', 'Cav3_1', 'pas')
}
Purkinje.ion_channels['Kv1_5'] = NMODLChannel('Kv1_5', conductance_density_parameter="gKur")
Purkinje.ion_channels['Cav3_1'] = NMODLChannel('Cav3_1', conductance_density_parameter="pcabar")
Purkinje.ion_channels['Cav2_1'] = NMODLChannel('Cav2_1', conductance_density_parameter="pcabar")
Purkinje.ion_channels['Cav3_3'] = NMODLChannel('Cav3_3', conductance_density_parameter="pcabar")
Purkinje.ion_channels['cdp5'] = NMODLChannel('cdp5', conductance_density_parameter="Nannuli")  # not an ion channel, not conductance_density
                                                                                               # need to rethink how we decide where to insert channels

# Load cm values
subsets_cm = np.genfromtxt("ModelViewParmSubset_cm.txt")
cm_distributions = [
    uniform(with_label("soma", "axonAIS", "axonAISK", "axonNOR", "axonNOR2", "axonNOR3"), 0.77),
    uniform(with_label("axonmyelin", "axonmyelin2", "axonmyelin3", "axonmyelin4"), 1.87e-11),
    uniform(with_label("axoncoll", "axoncoll2"), 1.0),
    uniform(with_label("b0s02[24]"), 8.58298 * 0.77/1.64)
]
for subset_id, value in subsets_cm:
    cm_distributions.append(
        uniform(with_label("dend_subset{}".format(int(subset_id))), value * 0.77/1.64)
    )

purkinje_cell = Purkinje(
    morphology=morph,
    cm=any(*cm_distributions, absence=0.77),
    Ra=122,
    ionic_species={
        "h": IonicSpecies("h", reversal_potential=-34.4),
        "na": IonicSpecies("na", reversal_potential=uniform(with_label("axonAIS"), 75, absence=60)),
        "k": IonicSpecies("k", reversal_potential=-88),
        "ca": IonicSpecies("ca", 
                           reversal_potential=uniform(
                               with_label("soma", "dend", "axonAIS", "axonAISK", "axoncoll", "axoncoll2"),
                               137.52625, absence=None),
                           internal_concentration= 5e-5, external_concentration=2.0)
    },
    **ion_channel_parameters
)

population = sim.Population(1, purkinje_cell)
