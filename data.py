delta_o_cm = {
    "[Ti(H2O)6]3+": 20300,
    "[V(H2O)6]2+": 12600,
    "[V(H2O)6]3+": 18900,
    "[Cr(H2O)6]2+": 13900,
    "[CrCl6]3-": 13000,
    "[Cr(H2O)6]3+": 17400,
    "[Cr(NH3)6]3+": 21500,
    "[Cr(CN)6]3-": 26600,
    "Cr(CO)6": 34150,
    "[MnCl6]4-": 7500,
    "[Mn(H2O)6]2+": 8500,
    "[MnCl6]3-": 20000,
    "[Mn(H2O)6]3+": 21000,
    "[Fe(H2O)6]2+": 10400,
    "[Fe(H2O)6]3+": 14300,
    "[Fe(CN)6]4-": 32800,
    "[Fe(CN)6]3-": 35000,
    "[CoF6]3-": 13000,
    "[Co(H2O)6]2+": 9300,
    "[Co(H2O)6]3+": 27000,
    "[Co(NH3)6]3+": 22900,
    "[Co(CN)6]3-": 34800,
    "[Ni(H2O)6]2+": 8500,
    "[Ni(NH3)6]2+": 10800,
    "[RhCl6]3-": 20400,
    "[Rh(H2O)6]3+": 27000,
    "[Rh(NH3)6]3+": 34000,
    "[Rh(CN)6]3-": 45500,
    "[IrCl6]3-": 25000,
    "[Ir(NH3)6]3+": 41000,
}

ligand_data = {
    "Br-": {"f": 0.72, "dent": 1, "LF": 0.72},
    "(SCN)-": {"f": 0.75, "dent": 1, "LF": 0.75},
    "Cl-": {"f": 0.78, "dent": 1, "LF": 0.78},
    "OPCl3": {"f": 0.82, "dent": 1, "LF": 0.82},
    "(N3)-": {"f": 0.83, "dent": 1, "LF": 0.83},
    "F-": {"f": 0.90, "dent": 1, "LF": 0.90},
    "(OSMe3)-": {"f": 0.91, "dent": 1, "LF": 0.91},
    "OCMe2": {"f": 0.92, "dent": 1, "LF": 0.92},
    "EtOH": {"f": 0.97, "dent": 1, "LF": 0.97},
    "Me2NCHO": {"f": 0.98, "dent": 1, "LF": 0.98},
    "(C2O4)2-": {"f": 0.99, "dent": 2, "LF": 0.99},            # oxalate, O,O‑bidentate
    "H2O": {"f": 1.00, "dent": 1, "LF": 1.00},
    "SC(NH2)2": {"f": 1.01, "dent": 1, "LF": 1.01},
    "NCS-": {"f": 1.02, "dent": 1, "LF": 1.02},
    "NCSe-": {"f": 1.03, "dent": 1, "LF": 1.03},
    "(NH2CH2CO2)-": {"f": 1.18, "dent": 2, "LF": 1.18},        # glycinate, N,O‑bidentate
    "CH3CN": {"f": 1.22, "dent": 1, "LF": 1.22},
    "C5H5": {"f": 1.23, "dent": 1, "LF": 1.23},                # η5‑Cp treated as one site
    "NH3": {"f": 1.25, "dent": 1, "LF": 1.33},
    "en": {"f": 1.28, "dent": 2, "LF": 1.40},                  # ethylenediamine, N,N‑bidentate
    "(SO3)2-": {"f": 1.20, "dent": 2, "LF": 1.20},             # sulfite, O,O‑bidentate
    "diars": {"f": 1.33, "dent": 2, "LF": 1.33},               # diarsine, As,As‑bidentate
    "bipy": {"f": 1.33, "dent": 2, "LF": 1.33},                # 2,2′‑bipyridine, N,N‑bidentate
    "(NO2)-": {"f": 1.40, "dent": 1, "LF": 1.40},
    "CN-": {"f": 1.70, "dent": 1, "LF": 1.51},
}

metal_g = {
    "Ti2+": 695,
    "V2+": 755,
    "V3+": 861,
    "Cr2+": 810,
    "Cr3+": 918,
    "Mn2+": 860,
    "Mn3+": 965,
    "Fe2+": 917,
    "Fe3+": 1015,
    "Co2+": 971,
    "Co3+": 1065,
    "Ni2+": 1030,
    "Ni3+": 1115,
}

# ---------------------------------------------
# Empirical 10 Dq values (cm⁻¹) for hexaaqua ions
# ---------------------------------------------
metal_DqH2O = {
    "Cr2+": 13900,  # Cr(H2O)6 2+
    "Cr3+": 17400,  # Cr(H2O)6 3+
    "V3+":  18900,  # V(H2O)6 3+
    "Co2+":  9300,  # Co(H2O)6 2+
    "Co3+": 18600,  # Co(H2O)6 3+
    "Ni2+":  8600,  # Ni(H2O)6 2+
}

geometry_coord = {
    "Octahedral": 6,
    "Tetrahedral": 4,
    "Square-Pyramidal": 5,
    "Trigonal-Bipyramidal": 5,
    "Square-Planar": 4,
    "Trigonal-Planar": 3,
    "Linear": 2,
}

# empirical geometry scaling factors (relative to octahedral Δ)
geometry_scale = {
    "Octahedral": 1.00,
    "Tetrahedral": 4/9,          # ≈0.44
    "Square-Pyramidal": 0.95,
    "Trigonal-Bipyramidal": 0.85,
    "Square-Planar": 1.35,       # typical for d8 metals
    "Trigonal-Planar": 0.70,
    "Linear": 0.30,
}

# ---------------------------------------------
# Electronic configuration and Tanabe–Sugano multipliers
# ---------------------------------------------
# For each metal entry in `metal_g` we record its d‑electron count and a
# representative spin state (HS = high‑spin, LS = low‑spin).  The second
# dictionary gives an empirical multiplier that converts Δ (10 Dq) into
# the energy of the *first* intense spin‑allowed transition, based on
# Tanabe–Sugano charts.  These factors are approximate but bring the
# predicted λmax into the correct visible region for common ions.

electron_config = {
    "Ti2+": ("d2", "HS"),
    "V2+":  ("d3", "HS"),
    "V3+":  ("d2", "HS"),
    "Cr2+": ("d4", "HS"),
    "Cr3+": ("d3", "HS"),
    "Mn2+": ("d5", "HS"),
    "Mn3+": ("d4", "HS"),
    "Fe2+": ("d6", "HS"),
    "Fe3+": ("d5", "LS"),   # Fe³⁺ often low‑spin with strong‑field ligands
    "Co2+": ("d7", "HS"),
    "Co3+": ("d6", "LS"),
    "Ni2+": ("d8", "HS"),
    "Ni3+": ("d7", "LS"),
}

# (dn, spin)  → multiplier applied to Δ to estimate lowest‑energy
# spin‑allowed absorption. Values rounded from the TS diagrams.
ts_multiplier = {
    ("d1", "HS"): 1.00,
    ("d2", "HS"): 1.20,
    ("d3", "HS"): 1.70,
    ("d4", "HS"): 1.80,
    ("d5", "HS"): 1.50,
    ("d6", "HS"): 0.80,
    ("d7", "HS"): 1.20,
    ("d8", "HS"): 1.70,
    ("d4", "LS"): 1.10,
    ("d5", "LS"): 2.20,
    ("d6", "LS"): 1.00,
    ("d7", "LS"): 1.20,
}

# ---------------------------------------------
# Tanabe–Sugano spin‑allowed band ratios
# Each list entry is E_band / (10 Dq) for that d‑electron config.
# ---------------------------------------------
ts_bands = {
    # Octahedral (Oh)
    ("d2", "HS", "Octahedral"): [1.00, 1.49],
    ("d3", "HS", "Octahedral"): [1.00, 1.50],
    ("d4", "HS", "Octahedral"): [1.00],
    ("d6", "LS", "Octahedral"): [1.00, 1.32],
    ("d7", "HS", "Octahedral"): [1.20, 1.50],
    ("d8", "HS", "Octahedral"): [1.00, 1.63, 2.62],
    # Tetrahedral (Td)
    ("d7", "HS", "Tetrahedral"): [4.60],
}

# Charge‑transfer band energies (cm⁻¹) for strong π‑acceptor ligands
ct_band = {
    "CN-": 24000,   # ~417 nm
}