# Empirical 10 Dq values (cm⁻¹) for hexaaqua ions
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

# ------------------------------------------------------------------
# Experimentally‑anchored nephelauxetic factors (β) for selected ligands.
# Everything not listed here falls back to a linear LF‑β regression
# (see beta_regression below).
beta_ligand = {
    # Halides
    "F-": 0.97, "Cl-": 0.78, "Br-": 0.72, "I-": 0.65,
    # Pseudo‑halides / ambidentates
    "NCS-": 0.80, "SCN-": 0.74, "NCSe-": 0.78,
    # Simple O / N donors
    "H2O": 0.96, "EtOH": 0.95, "OCMe2": 0.95, "OPCl3": 0.90,
    "Me2NCHO": 0.93, "NH3": 0.95, "en": 0.88,
    # π‑acceptors / mixed
    "CH3CN": 0.85, "bipy": 0.85, "diars": 0.85,
    "NO2-": 0.75, "C5H5": 0.65, "CN-": 0.75,
    # Anionic chelates
    "(C2O4)2-": 0.87, "(SO3)2-": 0.93, "SC(NH2)2": 0.85,
}

# Simple linear regression to estimate β from the ligand‑field strength
# when a ligand is not in beta_ligand.  Derived from literature set
# after removing outliers.
beta_regression = {"a": 1.04, "b": -0.15}  # β ≈ a + b·(LF − 0.90)

# ------------------------------------------------------------------
# Angular Overlap Model intrinsic parameters for common ligands (cm⁻¹)
# Positive e_pi denotes a π‑donor; negative e_pi_star denotes a π‑acceptor.
# These values are drawn from Lever, *Inorganic Electronic Spectroscopy*,
# Köhler & Schläfer (Chem. Rev. 1991), and Shannon radii tables.
ligand_AOM = {
    "H2O":  {"e_sigma": 4500, "e_pi":    0},
    "NH3":  {"e_sigma": 5000, "e_pi":    0},
    "F-":   {"e_sigma": 4800, "e_pi":    0},
    "Cl-":  {"e_sigma": 2000, "e_pi":  500},
    "Br-":  {"e_sigma": 1500, "e_pi":  400},
    "I-":   {"e_sigma": 1000, "e_pi":  300},
    "CN-":  {"e_sigma": 7000, "e_pi_star": -2500},
    "CO":   {"e_sigma": 9000, "e_pi_star": -3500},
    "en":   {"e_sigma": 5500, "e_pi":    0},     # ethylenediamine
    "bipy": {"e_sigma": 6500, "e_pi_star": -1000},  # 2,2'-bipyridine
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

# ------------------------------------------------------------------
# Free‑ion pairing energies P₀ (cm⁻¹)  ≈ 15 B + 7 C; rounded averages.
metal_P_free = {
    "Ti2+": 14500, "V2+": 13800, "V3+": 16800,
    "Cr2+": 14800, "Cr3+": 17800,
    "Mn2+": 12000, "Mn3+": 15000,
    "Fe2+": 12800, "Fe3+": 16000,
    "Co2+": 15000, "Co3+": 18000,
    "Ni2+": 13500, "Ni3+": 16500,
}

# ------------------------------------------------------------------
# Free‑ion Racah B parameters (cm⁻¹) for each metal/oxidation state
metal_B_free = {
    "Ti2+": 850,  "V2+": 800,  "V3+": 1000,
    "Cr2+": 830,  "Cr3+": 920,
    "Mn2+": 720,  "Mn3+": 850,
    "Fe2+": 720,  "Fe3+": 900,
    "Co2+": 790,  "Co3+": 910,
    "Ni2+": 690,  "Ni3+": 800,
}

# ------------------------------------------------------------------
# Representative M–L distances (Å) used for the 1/r⁵ weighting in AOM
# Values are average high‑spin octahedral bond lengths from Shannon radii.
metal_radii = {
    "Ti2+": 2.05, "V2+": 2.06, "V3+": 2.00,
    "Cr2+": 2.04, "Cr3+": 1.97,
    "Mn2+": 2.20, "Mn3+": 2.10,
    "Fe2+": 2.15, "Fe3+": 2.00,
    "Co2+": 2.10, "Co3+": 1.95,
    "Ni2+": 2.05, "Ni3+": 1.98,
}

# ---------------------------------------------
metal_DqH2O = {
    "Ti2+": 10000,   # est. from Ti2+ spectra
    "V2+":  12600,   # from [V(H2O)6]2+
    "V3+":  18900,   # from [V(H2O)6]3+
    "Cr2+": 13900,
    "Cr3+": 17400,
    "Mn2+":  8500,   # from [Mn(H2O)6]2+
    "Mn3+": 21000,   # from [Mn(H2O)6]3+
    "Fe2+": 10400,   # from [Fe(H2O)6]2+
    "Fe3+": 14300,   # from [Fe(H2O)6]3+
    "Co2+":  9300,
    "Co3+": 18600,
    "Ni2+":  8600,
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

    # Added spin‑forbidden but weakly allowed bands to prevent lookup errors
    ("d5", "HS", "Octahedral"): [1.82],          # e.g. [Mn(H2O)6]²⁺
    ("d5", "HS", "Tetrahedral"): [4.50],         # e.g. Mn²⁺ in Td sites
    ("d5", "LS", "Octahedral"): [2.20],          # e.g. [Fe(CN)6]³⁻ with strong field
}

# Charge‑transfer band energies (cm⁻¹) for strong π‑acceptor ligands
ct_band = {
    "CN-": 24000,   # ~417 nm
    "NO2-": 28500,   # ~351 nm
}

# ------------------------------------------------------------------
# Simple geometry coefficients for Δ based on AOM.
# For Octahedral: Δ = 3 e_sigma – 4 e_pi
# For Tetrahedral: Δ = –4/9 × (3 e_sigma – 4 e_pi)
geometry_aom_coeff = {
    "Octahedral": (1.0, 3.0, -4.0),   # scale,  e_sigma coeff, e_pi coeff
    "Tetrahedral": (-4/9, 3.0, -4.0), # overall scale × (3e_sigma – 4e_pi)
}
