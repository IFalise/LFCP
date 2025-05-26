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
    beta_ligand,       # NEW
    beta_regression,   # NEW
    metal_P_free,      # NEW
)
# ------------------------------------------------------------------
# Optional plotting backend (Plotly is lighter + Streamlit‑native)
try:
    import plotly.graph_objects as go
    _PLOT_AVAILABLE = True
except ModuleNotFoundError:
    _PLOT_AVAILABLE = False


def wavelength_nm(dq_cm: float) -> float:
    """Convert Δ₀ in cm⁻¹ to λmax in nm."""
    return 1e7 / dq_cm


def nm_to_rgb(nm: float) -> str:
    """Visual mapping of 380‑780 nm → Hex colour using a CIE‑based approximation."""
    if nm < 380 or nm > 780:
        return "#FFFFFF"

    # Intensity attenuation near the vision limits
    if 380 <= nm < 420:
        attenuation = 0.3 + 0.7 * (nm - 380) / 40
    elif 700 < nm <= 780:
        attenuation = 0.3 + 0.7 * (780 - nm) / 80
    else:
        attenuation = 1.0

    # Un‐gamma‑corrected RGB
    if   380 <= nm < 440:
        r, g, b = -(nm - 440) / 60, 0.0, 1.0
    elif 440 <= nm < 490:
        r, g, b = 0.0, (nm - 440) / 50, 1.0
    elif 490 <= nm < 510:
        r, g, b = 0.0, 1.0, -(nm - 510) / 20
    elif 510 <= nm < 580:
        r, g, b = (nm - 510) / 70, 1.0, 0.0
    elif 580 <= nm < 645:
        r, g, b = 1.0, -(nm - 645) / 65, 0.0
    else:  # 645‑700
        r, g, b = 1.0, 0.0, 0.0

    gamma = 0.8
    def _corr(c: float) -> int:
        return 0 if c <= 0 else int(round((attenuation * (c ** gamma)) * 255))

    return f"#{_corr(r):02X}{_corr(g):02X}{_corr(b):02X}"


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


# ------------------------------------------------------------------
def beta_from_lf(lf: float) -> float:
    """Fallback β from LF via linear regression."""
    a = beta_regression["a"]
    b = beta_regression["b"]
    β = a + b * (lf - 0.90)
    return max(0.65, min(1.05, β))

def average_beta(ligand_counts: dict[str, int]) -> float:
    """Donor‑site‑weighted average nephelauxetic factor β̄."""
    total_sites = sum(ligand_data[lig]["dent"] * n
                      for lig, n in ligand_counts.items())
    if total_sites == 0:
        return 0.92   # benign default
    β_sum = 0.0
    for lig, n in ligand_counts.items():
        sites = ligand_data[lig]["dent"] * n
        lf    = ligand_data[lig]["LF"]
        β     = beta_ligand.get(lig, beta_from_lf(lf))
        β_sum += β * sites
    return β_sum / total_sites

def predict_spin_state(metal: str,
                       dn: str,
                       geometry: str,
                       delta_cm: float,
                       β_bar: float) -> tuple[str, int]:
    """
    Decide HS/LS based on Δ vs. pairing energy P = β̄ P₀.
    Returns (spin, n_unpaired).
    """
    if geometry not in ("Octahedral", "Square-Planar"):
        return "HS", {"d4":4,"d5":5,"d6":4,"d7":3}.get(dn, 0)

    P0 = metal_P_free[metal]
    P  = β_bar * P0

    # low‑spin window only relevant for d4–d7
    if dn in ("d4","d5","d6","d7"):
        spin = "LS" if delta_cm > P else "HS"
    else:
        spin = "HS"

    unpaired_map = {
        ("d4","HS"):4, ("d4","LS"):2,
        ("d5","HS"):5, ("d5","LS"):1,
        ("d6","HS"):4, ("d6","LS"):0,
        ("d7","HS"):3, ("d7","LS"):1,
    }
    n_u = unpaired_map.get((dn, spin), 0)
    return spin, n_u


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
    absorbed_hex = nm_to_rgb(λ_nm)          # colour of the absorbed band
    solution_hex = complement_hex(absorbed_hex)  # complementary solution colour

    st.markdown(f"**Preset selected:** `{preset_choice}`")
    st.markdown(f"10 Dq = `{dq} cm⁻¹`  →  λ<sub>max</sub> ≈ `{λ_nm:.0f} nm`")
    # Display absorbed vs. complementary solution colour
    st.markdown(
        "<div style='display:flex;gap:20px;'>"
        f"<div style='text-align:center'>"
        f"<div style='width:120px;height:50px;border:1px solid #000;background:{absorbed_hex}'></div>"
        "<small>Absorption</small></div>"
        f"<div style='text-align:center'>"
        f"<div style='width:120px;height:50px;border:1px solid #000;background:{solution_hex}'></div>"
        "<small>Solution</small></div>"
        "</div>",
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
        # --- Experimental baseline vs. AOM ----------------------------------------
        lf_bar = average_lf(ligand_counts)

        if metal in metal_DqH2O:
            # Empirical route: scale the known hexaaqua splitting by LF
            ten_Dq_base = metal_DqH2O[metal] * lf_bar
        else:
            # Fallback to Angular Overlap Model if no reference value exists
            e_sigma, e_pi, e_pi_star = average_aom(ligand_counts)
            delta_oct = 3 * e_sigma - 4 * e_pi
            ten_Dq_base = abs(delta_oct)

        # Geometry conversion (apply only once, after reference scaling)
        if geometry == "Tetrahedral":
            ten_Dq_base *= 4 / 9
        elif geometry == "Square-Planar":
            ten_Dq_base *= 1.35
        else:
            ten_Dq_base *= geometry_scale.get(geometry, 1.0)

        # --- Nephelauxetic correction & spin state ---------------------
        β_bar  = average_beta(ligand_counts)
        dn, _  = electron_config[metal]        # ignore preset spin
        spin, n_unpaired = predict_spin_state(metal, dn, geometry, ten_Dq_base, β_bar)

        ratios = ts_bands.get((dn, spin, geometry)) \
            or ts_bands.get((dn, spin, "Octahedral")) \
            or [1.00]   # fallback to a single‑ratio list to avoid None

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

        def _show_spectrum(nms: list[float]) -> None:
            if not _PLOT_AVAILABLE:
                return  # plotly not available → skip spectrum
            fig = go.Figure()
            for nm in nms:
                fig.add_shape(
                    type="line",
                    x0=nm, x1=nm,
                    y0=0,  y1=1,
                    line=dict(color=nm_to_rgb(nm), width=3),
                )
            fig.update_layout(
                xaxis=dict(range=[380, 780], title="λ (nm)", tickmode="array",
                           tickvals=[400, 500, 600, 700]),
                yaxis=dict(visible=False),
                showlegend=False,
                height=180,
                margin=dict(l=0, r=0, t=30, b=0),
                title_text="Predicted absorption bands",
                title_font_size=12,
            )
            st.plotly_chart(fig, use_container_width=False)

        # Show spectrum only if matplotlib is present (see _PLOT_AVAILABLE flag)
        _show_spectrum(visible_nms)

        st.markdown(
            f"**Spin state:** {spin} &nbsp;|&nbsp; **Unpaired e⁻:** {n_unpaired} "
            f"&nbsp;|&nbsp; β̄ ≈ {β_bar:.2f}"
        )

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
