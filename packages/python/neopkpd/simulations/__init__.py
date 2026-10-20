"""
NeoPKPD Simulations Package

This package provides simulation functions for various PK and PD models.

All PK models now support:
- Infusion administration via 'duration' parameter in dose events
- Absorption lag time (ALAG) via 'alag' parameter
- Bioavailability (F) via 'bioavailability' parameter
"""

from .pk_onecomp import (
    simulate_pk_iv_bolus,
    simulate_pk_oral_first_order,
)

from .pk_twocomp import (
    simulate_pk_twocomp_iv_bolus,
    simulate_pk_twocomp_oral,
)

from .pk_threecomp import (
    simulate_pk_threecomp_iv_bolus,
)

from .pk_advanced import (
    simulate_pk_transit_absorption,
    simulate_pk_michaelis_menten,
)

from .pk_custom import (
    simulate_pk_tmdd_custom,
    simulate_pk_parallel_absorption,
    simulate_pk_enterohepatic_recirculation,
    simulate_pk_autoinduction,
)

from .pkpd import (
    # Basic effect models
    simulate_pkpd_direct_emax,
    simulate_pkpd_sigmoid_emax,
    simulate_pkpd_biophase_equilibration,
    # Indirect response models (IRM)
    simulate_pkpd_indirect_response,  # IRM-III (inhibition of Kout)
    simulate_pkpd_irm1,               # IRM-I (inhibition of Kin)
    simulate_pkpd_irm2,               # IRM-II (stimulation of Kin)
    simulate_pkpd_irm4,               # IRM-IV (stimulation of Kout)
    # Advanced PD models
    simulate_pkpd_transit_compartment,
    simulate_pkpd_disease_progression,
    # Tolerance models
    simulate_pkpd_tolerance_counter_regulation,
    simulate_pkpd_receptor_regulation,
)

from .population import (
    simulate_population_iv_bolus,
    simulate_population_oral,
)

from .sensitivity import (
    run_sensitivity,
)

from .gsa import (
    run_sobol_sensitivity,
    run_morris_sensitivity,
    SobolResult,
    SobolIndex,
    MorrisResult,
    MorrisIndex,
)


__all__ = [
    # One-compartment PK
    "simulate_pk_iv_bolus",
    "simulate_pk_oral_first_order",
    # Two-compartment PK
    "simulate_pk_twocomp_iv_bolus",
    "simulate_pk_twocomp_oral",
    # Three-compartment PK
    "simulate_pk_threecomp_iv_bolus",
    # Advanced PK
    "simulate_pk_transit_absorption",
    "simulate_pk_michaelis_menten",
    # Custom PK models
    "simulate_pk_tmdd_custom",
    "simulate_pk_parallel_absorption",
    "simulate_pk_enterohepatic_recirculation",
    "simulate_pk_autoinduction",
    # PKPD - Basic effect models
    "simulate_pkpd_direct_emax",
    "simulate_pkpd_sigmoid_emax",
    "simulate_pkpd_biophase_equilibration",
    # PKPD - Indirect response models (IRM)
    "simulate_pkpd_indirect_response",  # IRM-III
    "simulate_pkpd_irm1",               # Inhibition of Kin
    "simulate_pkpd_irm2",               # Stimulation of Kin
    "simulate_pkpd_irm4",               # Stimulation of Kout
    # PKPD - Advanced models
    "simulate_pkpd_transit_compartment",
    "simulate_pkpd_disease_progression",
    # PKPD - Tolerance models
    "simulate_pkpd_tolerance_counter_regulation",
    "simulate_pkpd_receptor_regulation",
    # Population
    "simulate_population_iv_bolus",
    "simulate_population_oral",
    # Sensitivity (Local)
    "run_sensitivity",
    # Global Sensitivity Analysis
    "run_sobol_sensitivity",
    "run_morris_sensitivity",
    "SobolResult",
    "SobolIndex",
    "MorrisResult",
    "MorrisIndex",
]
