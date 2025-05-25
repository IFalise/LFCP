import streamlit as st
from data import (
    delta_o_cm,
    metal_g,
    ligand_data,
    ligand_AOM,
    metal_B_free,
    metal_radii,
    geometry_coord,
    geometry_scale,
    geometry_aom_coeff,
    electron_config,
    ts_multiplier,
    metal_DqH2O,   # kept for presets
    ts_bands,
    ct_band,
)


def wavelength_nm(dq_cm: float) -> float:
    """Convert Δ₀ in cm⁻¹ to λmax in nm."""
    return 1e7 / dq_cm


def nm_to_rgb(nm: float) -> str:
    """Visual mapping of 380‑780 nm → Hex colour."""
    bands = [
        (380, 400, "#6D007E"),
        (400, 420, "#6700B4"),
        (420, 440, "#3700E6"),
        (440, 460, "#0046FF"),
        (460, 480, "#00A9FF"),
        (480, 500, "#00FFFF"),
        (500, 520, "#00FF00"),
        (520, 540, "#5DFF00"),
        (540, 560, "#A2FF00"),
        (560, 580, "#E1FF00"),
        (580, 600, "#FFDF00"),
        (600, 620, "#FF9B00"),
        (620, 640, "#FF4E00"),
        (640, 660, "#F80000"),
        (660, 680, "#DC0000"),
        (680, 700, "#BF0000"),
        (700, 780, "#A10000"),
    ]

    for lo, hi, hexcode in bands:
        if lo <= nm < hi:
            return hexcode
    return "#FFFFFF"


def complement_hex(hexcode: str) -> str:
    """Return the complementary (solution) colour for a given absorbed hex."""
    r = 255 - int(hexcode[1:3], 16)
    g = 255 - int(hexcode[3:5], 16)
    b = 255 - int(hexcode[5:7], 16)
    return f"#{r:02X}{g:02X}{b:02X}"


def mix_hex(*codes: str) -> str:
    """Return average of ≥1 hex colours."""
    r = sum(int(c[1:3], 16) for c in codes) // len(codes)
    g = sum(int(c[3:5], 16) for c in codes) // len(codes)
    b = sum(int(c[5:7], 16) for c in codes) // len(codes)
    return f"#{r:02X}{g:02X}{b:02X}"


def average_lf(ligand_counts: dict[str, int]) -> float:
    """Donor‑site‑weighted average ligand‑field strength LF."""
    site_sum = lf_sum = 0.0
    for lig, n_lig in ligand_counts.items():
        entry = ligand_data[lig]
        sites = n_lig * entry["dent"]
        site_sum += sites
        lf_sum += sites * entry["LF"]
    return lf_sum / site_sum if site_sum else 1.0


def average_aom(ligand_counts: dict[str, int]) -> tuple[float, float, float]:
    """
    Return donor‑site‑weighted averages of eσ, eπ and eπ* for the selected
    ligand set.  Missing ligands default to eσ=4000, eπ=0.
    """
    site_sum = es_sum = epi_sum = epi_star_sum = 0.0
    for lig, n_lig in ligand_counts.items():
        sites = ligand_data[lig]["dent"] * n_lig
        site_sum += sites
        entry = ligand_AOM.get(lig, {"e_sigma": 4000})
        es_sum += sites * entry.get("e_sigma", 0)
        epi_sum += sites * entry.get("e_pi", 0)
        epi_star_sum += sites * entry.get("e_pi_star", 0)
    if site_sum == 0:
        return 0, 0, 0
    return es_sum / site_sum, epi_sum / site_sum, epi_star_sum / site_sum


st.title("Ligand Field Colour Predictor")

preset_label = "— build custom complex —"
preset_choice = st.selectbox(
    "Choose a preset complex *or* build your own",
    [preset_label] + list(delta_o_cm.keys()),
)

# ── Preset complex ────────────────────────────────────────────────────────────
if preset_choice != preset_label:
    dq = delta_o_cm[preset_choice]
    λ_nm = wavelength_nm(dq)
    colour = nm_to_rgb(λ_nm)

    st.markdown(f"**Preset selected:** `{preset_choice}`")
    st.markdown(f"10 Dq = `{dq} cm⁻¹`  →  λ<sub>max</sub> ≈ `{λ_nm:.0f} nm`")
    st.markdown(
        f"<div style='width:120px;height:50px;border:1px solid #000;"
        f"background:{colour}'></div>",
        unsafe_allow_html=True,
    )

# ── Custom complex builder ────────────────────────────────────────────────────
else:
    st.subheader("Custom complex builder")

    # metal and geometry
    metal = st.selectbox("Metal / oxidation state", list(metal_g.keys()))
    geometry = st.radio("Geometry", list(geometry_coord.keys()), horizontal=True)
    max_sites = geometry_coord[geometry]

    # ligand table
    st.markdown("#### Ligand set")
    n_rows = st.number_input("Number of different ligands", 1, 6, 1, step=1)
    ligand_counts: dict[str, int] = {}
    total_sites = 0

    for i in range(int(n_rows)):
        col1, col2 = st.columns([3, 1])
        with col1:
            lig = st.selectbox(
                f"Ligand {i + 1}", list(ligand_data.keys()), key=f"lig_{i}"
            )
        with col2:
            cnt = st.number_input("count", 0, 6, step=1, key=f"cnt_{i}")
        ligand_counts[lig] = int(cnt)
        total_sites += ligand_data[lig]["dent"] * int(cnt)

    # live validation
    if total_sites > max_sites:
        st.error(
            f"Total donor sites = {total_sites} exceeds "
            f"{max_sites} for {geometry}."
        )
    elif total_sites < max_sites:
        st.warning(f"{max_sites - total_sites} open coordination sites remain.")
    else:
        st.success("Coordination number satisfied.")

    # submit button
    can_submit = (total_sites == max_sites) and any(ligand_counts.values())
    submit = st.button("Compute colour", disabled=not can_submit)

    # calculation
    if submit:
        # --- Angular Overlap Model ----------------------------------------
        e_sigma, e_pi, e_pi_star = average_aom(ligand_counts)

        # Base octahedral 10 Dq from AOM: ΔOct = 3 eσ – 4 eπ  (LibreTexts Eq {1})
        delta_oct = 3 * e_sigma - 4 * e_pi

        # Geometry conversion
        if geometry == "Octahedral":
            ten_Dq_base = delta_oct
        elif geometry == "Tetrahedral":
            ten_Dq_base = (-4 / 9) * delta_oct
        else:
            # fallback to existing empirical scale for less common geometries
            ten_Dq_base = delta_oct * geometry_scale.get(geometry, 1.0)

        ten_Dq_base = abs(ten_Dq_base)  # magnitude only for band energies

        dn, spin = electron_config[metal]
        ratios = ts_bands.get((dn, spin, geometry), ts_bands.get((dn, spin, "Octahedral")))

        visible_nms: list[float] = []
        for r in ratios:
            test_cm = r * ten_Dq_base
            nm = 1e7 / test_cm
            if 380 <= nm <= 780:
                visible_nms.append(nm)

        # fall back to first ratio if nothing visible
        if not visible_nms:
            visible_nms.append(1e7 / (ratios[0] * ten_Dq_base))

        # Charge‑transfer fallback (use intrinsic eπ* if available)
        if (len(visible_nms) == 1) and not (380 <= visible_nms[0] <= 780):
            if e_pi_star != 0:
                ct_cm = abs(-e_pi_star) + abs(delta_oct)
                ct_nm = 1e7 / ct_cm
                if 380 <= ct_nm <= 780:
                    visible_nms = [ct_nm]
            else:
                for lig in ligand_counts:
                    if lig in ct_band:
                        ct_nm = 1e7 / ct_band[lig]
                        if 380 <= ct_nm <= 780:
                            visible_nms = [ct_nm]
                            break

        # Colour swatch: mix complements of all visible bands
        absorbed_hexes = [nm_to_rgb(nm) for nm in visible_nms]
        colour = mix_hex(*(complement_hex(hx) for hx in absorbed_hexes))

        # Report the *lowest‑energy* visible band
        lam_nm = min(visible_nms)
        band_cm = 1e7 / lam_nm

        st.markdown(f"Metal: {metal} &nbsp; | &nbsp; Geometry: {geometry}")
        st.markdown(
            f"10 Dq ≈ `{band_cm:.0f} cm⁻¹`  →  λ<sub>max</sub> ≈ `{lam_nm:.0f} nm`"
        )
        st.markdown(
            f"<div style='width:120px;height:50px;border:1px solid #000;"
            f"background:{colour}'></div>",
            unsafe_allow_html=True,
        )
