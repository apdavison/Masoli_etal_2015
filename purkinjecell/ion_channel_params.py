"""

"""

from __future__ import division
from math import exp
from pyNN.morphology import with_label, uniform, by_diameter, any


ion_channel_parameters = {
    "Leak": {
        "e": uniform(with_label("soma", "dend", "axonAIS", "axonAISK", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"),
                     -63),
        "gmax": any(
            uniform(with_label("soma"), 1.1e-3),
            uniform(with_label("axonAIS", "axonAISK", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 0.0003),
            uniform(with_label("b0s02[24]"), 1.74451e-4 / 2), 
            uniform(with_label("dend_subset10"), 9.81576e-5 / 2),
            uniform(with_label("dend_subset11"), 9.93235e-5 / 2),
            uniform(with_label("dend_subset14"), 1.03622e-4 / 2),
            uniform(with_label("dend_subset7"), 9.44842e-5 / 2),
            uniform(with_label("dend_subset48"), 3.07143e-4 / 2),
            uniform(with_label("dend_subset6"), 9.32852e-5 / 2),
            uniform(with_label("dend_subset44"), 2.28496e-4 / 2),
            uniform(with_label("dend_subset8"), 9.57322e-5 / 2),
            uniform(with_label("dend_subset9"), 9.68117e-5 / 2),
            uniform(with_label("dend_subset13"), 1.02042e-4 / 2),
            uniform(with_label("dend_subset26"), 1.27897e-4 / 2),
            uniform(with_label("dend_subset5"), 9.23213e-5 / 2),
            uniform(with_label("dend_subset30"), 1.39992e-4 / 2),
            uniform(with_label("dend_subset45"), 2.3946e-4 / 2),
            uniform(with_label("dend_subset24"), 1.22388e-4 / 2),
            uniform(with_label("dend_subset21"), 1.1597e-4 / 2),
            uniform(with_label("dend_subset22"), 1.17874e-4 / 2),
            uniform(with_label("dend_subset15"), 1.04995e-4 / 2),
            uniform(with_label("dend_subset47"), 2.68529e-4 / 2),
            uniform(with_label("dend_subset34"), 1.55635e-4 / 2),
            uniform(with_label("dend_subset46"), 2.54361e-4 / 2),
            uniform(with_label("dend_subset17"), 1.08519e-4 / 2),
            uniform(with_label("dend_subset38"), 1.76656e-4 / 2),
            uniform(with_label("dend_subset32"), 1.47279e-4 / 2),
            uniform(with_label("dend_subset36"), 1.65314e-4 / 2),
            uniform(with_label("dend_subset25"), 1.2506e-4 / 2),
            uniform(with_label("dend_subset42"), 2.06402e-4 / 2),
            uniform(with_label("dend_subset41"), 1.98606e-4 / 2),
            uniform(with_label("dend_subset19"), 1.12068e-4 / 2),
            uniform(with_label("dend_subset12"), 1.00779e-4 / 2),
            uniform(with_label("dend_subset29"), 1.36397e-4 / 2),
            uniform(with_label("dend_subset33"), 1.50931e-4 / 2),
            uniform(with_label("dend_subset31"), 1.43185e-4 / 2),
            uniform(with_label("dend_subset37"), 1.71268e-4 / 2),
            uniform(with_label("dend_subset43"), 2.16786e-4 / 2),
            uniform(with_label("dend_subset40"), 1.9013e-4 / 2),
            uniform(with_label("dend_subset27"), 1.30398e-4 / 2),
            uniform(with_label("dend_subset23"), 1.20278e-4 / 2),
            uniform(with_label("dend_subset39"), 1.83704e-4 / 2),
            uniform(with_label("dend_subset18"), 1.10092e-4 / 2),
            uniform(with_label("dend_subset16"), 1.06714e-4 / 2),
            uniform(with_label("dend_subset35"), 1.60731e-4 / 2),
            uniform(with_label("dend_subset28"), 1.33581e-4 / 2),
            uniform(with_label("dend_subset20"), 1.13795e-4 / 2),
            uniform(with_label("dend_subset87"), 3.33333e-5 / 2),
            uniform(with_label("dend"), 3.33333e-5 / 2)
        )
    },
    "Cav3_1": {
        "pcabar": any(
            uniform(with_label("soma"), 7e-6),
            by_diameter(with_label("dend"), lambda d: (3.5 <= d <=12) and 5e-6 or 0.0),
            uniform(with_label("axonAIS"), 8.2e-6),
            uniform(with_label("axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 1e-5)
        )
    },
    "Cav2_1": {
        "pcabar": any(
            uniform(with_label("soma", "axonAIS", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 2.2e-4),
            uniform(with_label("dend"), 1e-3)
        )
    }, 
    "HCN1": {
        "gbar": any(
            uniform(with_label("soma"), 0.0004),
            uniform(with_label("dend"), 0.000004)
        )
    },
    "Nav1_6": {
        "gbar": any(
            uniform(with_label("soma"), 0.214),
            by_diameter(with_label("dend"), lambda d: (8 <= d <=12) and 0.016 or 0.0),
            uniform(with_label("axonAIS"), 0.50),
            uniform(with_label("axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 0.03)
        )
    },
    "Kv3_4": {
        "gkbar": any(
            uniform(with_label("soma"), 0.05),
            uniform(with_label("axonAIS"), 0.01),
            uniform(with_label("axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 0.02)
        )
    },
    "Kv1_1": {
        "gbar": any(
            uniform(with_label("soma"), 0.002),
            uniform(with_label("dend"), 0.0012),
            uniform(with_label("axonAISK"), 0.01)
        )
    }, 
    "Cav3_2": {
        "gcabar": any(
            uniform(with_label("soma"), 0.0008),
            by_diameter(with_label("dend"), lambda d: (3.5 <= d <=12) and 0.0012 or 0.0)
        )
    }, 
    "Kca3_1": {
        "gkbar": any(
            uniform(with_label("soma"), 0.01),
            by_diameter(with_label("dend"), lambda d: (3.5 <= d <=12) and 0.002 or 0.0)
        )
    },
    "Cav3_3": {
        "pcabar": uniform(with_label("soma", "dend"), 0.0001)
    }, 
    "Kir2_3": {
        "gkbar": any(
            uniform(with_label("soma"), 0.00003),
            by_diameter(with_label("dend"), lambda d: (3.5 <= d <=12) and 0.00001 or 0.0)
        )
    },
    "Kca1_1": {
        "gbar": any(
            uniform(with_label("soma"), 0.01),
            uniform(with_label("dend"), 3.5e-2)
        )
    },
    "Kca2_2": {
        "gkbar": any(
            uniform(with_label("soma"), 1e-3),
            by_diameter(with_label("dend"), lambda d: (3.5 <= d <=12) and 1e-3 or 0.0)
        )
    },
    "Kv4_3": {
        "gkbar": uniform(with_label("dend"), 0.001)
    },
    "Kv1_5": {
        "gKur": uniform(with_label("dend"), 0.13195e-3)
    },
    "Kv3_3": {
        "gbar": uniform(with_label("dend"), 0.01)
    },
    "cdp5": {
        "TotalPump": any(
            uniform(with_label("soma", "axonAIS"), 5e-8),
            uniform(with_label("dend"), 2e-8),
            uniform(with_label("axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 5e-7)
        ),
        "Nannuli": by_diameter(with_label("soma", "dend", "axonAIS", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"),
                               lambda d: 0.326 + (1.94 * d) + (0.289 * d**2) - (3.33e-2 * d**3) + (1.55e-3 * d**4) - (2.55e-5 * d**5)),
        "Buffnull2": by_diameter(with_label("soma", "dend", "axonAIS", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"),
                                 lambda d: 64.2 - 57.3 * exp(-d/1.4)),
        "rf3": by_diameter(with_label("soma", "dend", "axonAIS", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"),
                                      lambda d: 0.162 - 0.106 * exp(-d/2.29)),
        "rf4": any(
            by_diameter(with_label("soma"), 
                        lambda d: 0.000267 + 0.0167 * exp(-d/0.722) + 0.0028 * exp(-d/4)),
            uniform(with_label("axonAIS", "axonNOR", "axonNOR2", "axonNOR3", "axoncoll", "axoncoll2"), 0.003),
            by_diameter(with_label("dend"), 
                        lambda d: (d >= 2) and 0.000267 + 0.0167 * exp(-d/0.722) + 0.0028 * exp(-d/4) or 0.003)
        )
    },
    "pas": {
        "e": uniform(with_label("axonmyelin", "axonmyelin2", "axonmyelin3", "axonmyelin4"), -63),
        "g": uniform(with_label("axonmyelin", "axonmyelin2", "axonmyelin3", "axonmyelin4"), 5.6e-9)
    }
}