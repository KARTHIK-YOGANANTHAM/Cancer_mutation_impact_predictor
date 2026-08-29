"""
MUTATION SEVERITY INDEX (MSI) ANALYZER - ENTERPRISE EDITION v3.0
================================================================================
A high-precision computational bioinformatics platform for scoring the structural
and evolutionary severity of non-synonymous single amino acid substitutions.

Refined Visualizations, Expanded Domain Annotations, and Grantham Metric Engine.
================================================================================
"""

import bisect
import json
import math
import streamlit as st
import plotly.graph_objects as go

# ------------------------------------------------------------------------------
# 1. REFERENCE DATASETS & BIOINFORMATIC CONSTANTS
# ------------------------------------------------------------------------------

AA_PROPERTIES = {
    "A": {"name": "Alanine",       "hydrophobicity":  1.8, "charge":  0, "volume":  88.6, "category": "Nonpolar aliphatic"},
    "R": {"name": "Arginine",      "hydrophobicity": -4.5, "charge":  1, "volume": 173.4, "category": "Basic (positive)"},
    "N": {"name": "Asparagine",    "hydrophobicity": -3.5, "charge":  0, "volume": 114.1, "category": "Polar uncharged"},
    "D": {"name": "Aspartate",     "hydrophobicity": -3.5, "charge": -1, "volume": 111.1, "category": "Acidic (negative)"},
    "C": {"name": "Cysteine",      "hydrophobicity":  2.5, "charge":  0, "volume": 108.5, "category": "Polar uncharged"},
    "Q": {"name": "Glutamine",     "hydrophobicity": -3.5, "charge":  0, "volume": 143.8, "category": "Polar uncharged"},
    "E": {"name": "Glutamate",     "hydrophobicity": -3.5, "charge": -1, "volume": 138.4, "category": "Acidic (negative)"},
    "G": {"name": "Glycine",       "hydrophobicity": -0.4, "charge":  0, "volume":  60.1, "category": "Special (flexible)"},
    "H": {"name": "Histidine",     "hydrophobicity": -3.2, "charge":  1, "volume": 153.2, "category": "Basic (positive)"},
    "I": {"name": "Isoleucine",    "hydrophobicity":  4.5, "charge":  0, "volume": 166.7, "category": "Nonpolar aliphatic"},
    "L": {"name": "Leucine",       "hydrophobicity":  3.8, "charge":  0, "volume": 166.7, "category": "Nonpolar aliphatic"},
    "K": {"name": "Lysine",        "hydrophobicity": -3.9, "charge":  1, "volume": 168.6, "category": "Basic (positive)"},
    "M": {"name": "Methionine",    "hydrophobicity":  1.9, "charge":  0, "volume": 162.9, "category": "Nonpolar (sulfur)"},
    "F": {"name": "Phenylalanine", "hydrophobicity":  2.8, "charge":  0, "volume": 189.9, "category": "Aromatic"},
    "P": {"name": "Proline",       "hydrophobicity": -1.6, "charge":  0, "volume": 112.7, "category": "Special (rigid)"},
    "S": {"name": "Serine",        "hydrophobicity": -0.8, "charge":  0, "volume":  89.0, "category": "Polar uncharged"},
    "T": {"name": "Threonine",     "hydrophobicity": -0.7, "charge":  0, "volume": 116.1, "category": "Polar uncharged"},
    "W": {"name": "Tryptophan",    "hydrophobicity": -0.9, "charge":  0, "volume": 227.8, "category": "Aromatic"},
    "Y": {"name": "Tyrosine",      "hydrophobicity": -1.3, "charge":  0, "volume": 193.6, "category": "Aromatic (polar)"},
    "V": {"name": "Valine",        "hydrophobicity":  4.2, "charge":  0, "volume": 140.0, "category": "Nonpolar aliphatic"},
}

# Grantham Parameters: composition (c), polarity (p), molecular volume (v)
GRANTHAM_PROPS = {
    "A": {"c": 0.00, "p": 8.1,  "v": 31.0},  "R": {"c": 0.65, "p": 10.5, "v": 124.0},
    "N": {"c": 1.33, "p": 11.6, "v": 56.0},  "D": {"c": 1.33, "p": 13.0, "v": 54.0},
    "C": {"c": 0.00, "p": 5.5,  "v": 55.0},  "Q": {"c": 1.00, "p": 10.5, "v": 85.0},
    "E": {"c": 1.00, "p": 12.3, "v": 83.0},  "G": {"c": 0.00, "p": 9.0,  "v": 3.0},
    "H": {"c": 0.50, "p": 10.4, "v": 96.0},  "I": {"c": 0.00, "p": 5.2,  "v": 111.0},
    "L": {"c": 0.00, "p": 4.9,  "v": 111.0}, "K": {"c": 0.50, "p": 11.3, "v": 119.0},
    "M": {"c": 0.00, "p": 5.7,  "v": 105.0}, "F": {"c": 0.00, "p": 5.2,  "v": 132.0},
    "P": {"c": 0.00, "p": 8.0,  "v": 32.5},  "S": {"c": 1.00, "p": 9.2,  "v": 32.0},
    "T": {"c": 0.50, "p": 8.6,  "v": 61.0},  "W": {"c": 0.25, "p": 5.4,  "v": 170.0},
    "Y": {"c": 0.25, "p": 6.2,  "v": 136.0}, "V": {"c": 0.00, "p": 5.9,  "v": 84.0},
}

_AA_ORDER = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
_BLOSUM62_ROWS = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4],
]
BLOSUM62 = {(_AA_ORDER[i], _AA_ORDER[j]): _BLOSUM62_ROWS[i][j] for i in range(20) for j in range(20)}

GENE_INFO = {
    "KRAS": {
        "full_name": "Kirsten Rat Sarcoma Viral Oncogene Homolog",
        "function": "Small GTPase cycling between active GTP-bound and inactive GDP-bound states to mediate MAPK/ERK signaling.",
        "protein_length": 188,
        "domains": [
            {"name": "P-loop / G1 box", "start": 10, "end": 17},
            {"name": "Switch I", "start": 25, "end": 40},
            {"name": "Switch II", "start": 57, "end": 76},
            {"name": "Hypervariable Region", "start": 166, "end": 188},
        ],
    },
    "TP53": {
        "full_name": "Tumor Protein P53",
        "function": "Nuclear transcription factor orchestrating DNA repair, cell-cycle arrest, and apoptosis upon genomic stress.",
        "protein_length": 393,
        "domains": [
            {"name": "Transactivation Domain", "start": 1, "end": 42},
            {"name": "Proline-Rich Domain", "start": 63, "end": 97},
            {"name": "DNA-Binding Domain", "start": 102, "end": 292},
            {"name": "Tetramerization Domain", "start": 323, "end": 356},
        ],
    },
    "EGFR": {
        "full_name": "Epidermal Growth Factor Receptor",
        "function": "Receptor tyrosine kinase initiating RAS-RAF-MEK-ERK and PI3K-AKT cell growth cascades upon ligand binding.",
        "protein_length": 1210,
        "domains": [
            {"name": "Extracellular Ligand-Binding", "start": 1, "end": 620},
            {"name": "Transmembrane Domain", "start": 621, "end": 644},
            {"name": "Tyrosine Kinase Domain", "start": 712, "end": 979},
            {"name": "C-Terminal Regulatory Tail", "start": 980, "end": 1210},
        ],
    },
    "BRAF": {
        "full_name": "B-Raf Proto-Oncogene, Serine/Threonine Kinase",
        "function": "Downstream effector of RAS in the canonical RAF-MEK-ERK signaling pathway.",
        "protein_length": 766,
        "domains": [
            {"name": "RAS-Binding Domain (RBD)", "start": 155, "end": 227},
            {"name": "Cysteine-Rich Domain (CRD)", "start": 234, "end": 280},
            {"name": "Protein Kinase Domain", "start": 457, "end": 717},
        ],
    },
    "PIK3CA": {
        "full_name": "Phosphatidylinositol-4,5-Bisphosphate 3-Kinase Catalytic Subunit Alpha",
        "function": "Catalytic subunit p110alpha converting PIP2 to PIP3 in lipid signaling.",
        "protein_length": 1068,
        "domains": [
            {"name": "Adaptor-Binding Domain", "start": 32, "end": 108},
            {"name": "Ras-Binding Domain", "start": 180, "end": 290},
            {"name": "C2 Domain", "start": 350, "end": 485},
            {"name": "Helical Domain", "start": 520, "end": 696},
            {"name": "Kinase Domain", "start": 795, "end": 1080},
        ],
    },
    "PTEN": {
        "full_name": "Phosphatase and Tensin Homolog",
        "function": "Tumor suppressor lipid phosphatase counteracting PI3K signaling by converting PIP3 to PIP2.",
        "protein_length": 403,
        "domains": [
            {"name": "Phosphatase Catalytic Domain", "start": 14, "end": 185},
            {"name": "C2 Domain", "start": 186, "end": 351},
            {"name": "C-Terminal Tail", "start": 352, "end": 403},
        ],
    },
}

PRESET_MUTATIONS = {
    "KRAS G12D | Pancreatic & Colorectal Driver": {
        "gene": "KRAS", "position": 12, "wt": "G", "mut": "D",
        "note": "A major oncogenic codon substitution disrupting GTP hydrolysis and constitutive signaling.",
    },
    "TP53 R175H | DNA-Binding Domain Hotspot": {
        "gene": "TP53", "position": 175, "wt": "R", "mut": "H",
        "note": "Structural hotspot substitution inducing core domain destabilization and loss of DNA binding.",
    },
    "EGFR L858R | Lung Adenocarcinoma Driver": {
        "gene": "EGFR", "position": 858, "wt": "L", "mut": "R",
        "note": "Sensitizing mutation in the activation loop locking the kinase in a active state.",
    },
    "BRAF V600E | Melanoma Kinase Activation": {
        "gene": "BRAF", "position": 600, "wt": "V", "mut": "E",
        "note": "Mimics activation loop phosphorylation, causing ligand-independent kinase activity.",
    },
    "PIK3CA E545K | Helical Domain Hotspot": {
        "gene": "PIK3CA", "position": 545, "wt": "E", "mut": "K",
        "note": "Relieves inhibitory interactions with the p85 regulatory subunit.",
    },
    "PTEN R130Q | Catalytic Phosphatase Loss": {
        "gene": "PTEN", "position": 130, "wt": "R", "mut": "Q",
        "note": "A primary active-site catalytic mutation ablating lipid phosphatase activity.",
    },
}

COLOR_BENIGN = "#10B981"      # Emerald
COLOR_MODERATE = "#F59E0B"    # Amber
COLOR_PATHOGENIC = "#EF4444"  # Crimson
COLOR_PRIMARY = "#06B6D4"     # Cyan

# ------------------------------------------------------------------------------
# 2. COMPUTATIONAL & METRIC LOGIC
# ------------------------------------------------------------------------------

def compute_grantham_distance(wt_aa: str, mut_aa: str) -> float:
    """Calculate the Grantham amino acid difference distance (Grantham, 1974)."""
    p1, p2 = GRANTHAM_PROPS[wt_aa], GRANTHAM_PROPS[mut_aa]
    dc = p1["c"] - p2["c"]
    dp = p1["p"] - p2["p"]
    dv = p1["v"] - p2["v"]
    # Formula constants: alpha=1.833, beta=0.1018, gamma=0.000399 (scaled factor)
    dist = math.sqrt(1.833 * (dc**2) + 0.1018 * (dp**2) + 0.000399 * (dv**2)) * 10.0
    return round(dist, 1)


def compute_mutation_score(wt_aa: str, mut_aa: str) -> dict:
    """Compute composite Mutation Severity Index (MSI) and sub-components."""
    wt, mut = AA_PROPERTIES[wt_aa], AA_PROPERTIES[mut_aa]

    blosum_score = BLOSUM62[(wt_aa, mut_aa)]
    blosum_penalty = max(0.0, min(1.0, (11 - blosum_score) / 15.0))

    hydro_delta = abs(mut["hydrophobicity"] - wt["hydrophobicity"])
    hydro_norm = min(1.0, hydro_delta / 9.0)

    charge_delta = abs(mut["charge"] - wt["charge"])
    charge_norm = min(1.0, charge_delta / 2.0)

    vol_delta = abs(mut["volume"] - wt["volume"])
    vol_norm = min(1.0, vol_delta / 167.7)

    grantham_dist = compute_grantham_distance(wt_aa, mut_aa)
    grantham_norm = min(1.0, grantham_dist / 215.0)

    # Weighted composite severity index (0 - 100)
    msi = 100.0 * (0.35 * blosum_penalty + 0.25 * hydro_norm + 0.20 * charge_norm + 0.10 * vol_norm + 0.10 * grantham_norm)

    return {
        "msi": round(msi, 1),
        "blosum_score": blosum_score,
        "blosum_penalty": round(blosum_penalty, 3),
        "hydro_delta": round(hydro_delta, 2),
        "hydro_norm": round(hydro_norm, 3),
        "charge_delta": charge_delta,
        "charge_norm": round(charge_norm, 3),
        "vol_delta": round(vol_delta, 1),
        "vol_norm": round(vol_norm, 3),
        "grantham_dist": grantham_dist,
        "grantham_norm": round(grantham_norm, 3),
    }


def get_verdict(msi: float):
    if msi < 33.0:
        return "Tolerated / Likely Benign", COLOR_BENIGN, "benign"
    elif msi < 66.0:
        return "Uncertain / Moderate Disruption", COLOR_MODERATE, "moderate"
    else:
        return "High Impact / Severe Disruption", COLOR_PATHOGENIC, "pathogenic"


def classify_substitution(wt_aa: str, mut_aa: str):
    wt_cat = AA_PROPERTIES[wt_aa]["category"]
    mut_cat = AA_PROPERTIES[mut_aa]["category"]
    return wt_cat, mut_cat, wt_cat == mut_cat


def get_domain_context(gene: str, position: int):
    info = GENE_INFO.get(gene)
    if not info:
        return "Unannotated / Custom Gene"
    for d in info["domains"]:
        if d["start"] <= position <= d["end"]:
            return d["name"]
    return "Inter-domain Region"


@st.cache_data
def get_all_substitution_scores():
    scores = []
    for wt in AA_PROPERTIES:
        for mut in AA_PROPERTIES:
            if wt != mut:
                scores.append(compute_mutation_score(wt, mut)["msi"])
    return sorted(scores)


def get_percentile(msi: float, all_scores: list) -> float:
    idx = bisect.bisect_right(all_scores, msi)
    return round(100.0 * idx / len(all_scores), 1)


def generate_report_text(record: dict) -> str:
    scores = record["scores"]
    gene, wt_aa, mut_aa, position = record["gene"], record["wt"], record["mut"], record["position"]
    verdict_text, _, _ = get_verdict(scores["msi"])
    wt_cat, mut_cat, same_cat = classify_substitution(wt_aa, mut_aa)
    domain = get_domain_context(gene, position)
    all_scores = get_all_substitution_scores()
    pct = get_percentile(scores["msi"], all_scores)

    return f"""================================================================================
CANCER MUTATION IMPACT PREDICTOR - ANALYSIS REPORT
================================================================================
Target Gene       : {gene}
Substitution      : p.{AA_PROPERTIES[wt_aa]['name']} {position} {AA_PROPERTIES[mut_aa]['name']} ({wt_aa}{position}{mut_aa})
Domain Context    : {domain}

COMPOSITE SEVERITY INDEX
--------------------------------------------------------------------------------
MSI Score         : {scores['msi']} / 100
Clinical Verdict  : {verdict_text}
Severity Rank     : Upper {100 - pct}% (Exceeds {pct}% of all 380 non-synonymous substitutions)

PHYSICOCHEMICAL & EVOLUTIONARY BREAKDOWN
--------------------------------------------------------------------------------
BLOSUM62 Matrix Score  : {scores['blosum_score']} (Penalty Norm: {scores['blosum_penalty']})
Grantham Distance      : {scores['grantham_dist']} / 215.0
Hydrophobicity Shift   : {scores['hydro_delta']} (KD Scale)
Net Charge Shift       : {scores['charge_delta']} e
Residue Volume Shift   : {scores['vol_delta']} A^3
Structural Class Change: {wt_cat} -> {mut_cat} ({'Same Category' if same_cat else 'Class Transition'})

SCORING FORMULA METRICS
--------------------------------------------------------------------------------
MSI = 100 * [0.35(BLOSUM_penalty) + 0.25(Hydro_norm) + 0.20(Charge_norm) 
             + 0.10(Vol_norm) + 0.10(Grantham_norm)]

================================================================================
NOTE: Computational heuristic report generated for research analysis.
================================================================================
"""

def generate_report_json(record: dict) -> str:
    return json.dumps(record, indent=2)

# ------------------------------------------------------------------------------
# 3. HIGH-PERFORMANCE PLOTLY VISUALIZATIONS
# ------------------------------------------------------------------------------

_FONT = "Inter, -apple-system, BlinkMacSystemFont, sans-serif"
_MONO = "JetBrains Mono, monospace"
_DARK_BG = "#0B0F17"
_SURFACE_BG = "#111827"
_BORDER_COLOR = "#1F2937"
_TEXT_COLOR = "#F3F4F6"
_TEXT_MUTED = "#9CA3AF"


def _base_layout(fig, height=360, **extra):
    """Safely merge default layout kwargs with user extra args without key duplication."""
    layout_kwargs = dict(
        height=height,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font=dict(color=_TEXT_COLOR, family=_FONT),
        margin=dict(l=40, r=30, t=40, b=40),
    )
    layout_kwargs.update(extra)
    fig.update_layout(**layout_kwargs)
    return fig


def make_gauge_chart(msi: float, color: str) -> go.Figure:
    fig = go.Figure(go.Indicator(
        mode="gauge+number",
        value=msi,
        number={"suffix": " / 100", "font": {"size": 38, "family": _MONO, "color": color}},
        gauge={
            "axis": {"range": [0, 100], "tickcolor": _TEXT_MUTED, "tickfont": {"color": _TEXT_MUTED, "size": 11}},
            "bar": {"color": color, "thickness": 0.25},
            "bgcolor": "rgba(0,0,0,0)",
            "borderwidth": 0,
            "steps": [
                {"range": [0, 33], "color": "rgba(16, 185, 129, 0.1)"},
                {"range": [33, 66], "color": "rgba(245, 158, 11, 0.1)"},
                {"range": [66, 100], "color": "rgba(239, 68, 68, 0.1)"},
            ],
            "threshold": {"line": {"color": color, "width": 3}, "thickness": 0.8, "value": msi},
        },
    ))
    return _base_layout(fig, height=260, margin=dict(l=30, r=30, t=20, b=10))


def make_radar_chart(entries: list) -> go.Figure:
    categories = ["BLOSUM62 Penalty", "Hydrophobicity Shift", "Charge Change", "Volume Change", "Grantham Distance"]
    categories_closed = categories + [categories[0]]

    fig = go.Figure()
    for e in entries:
        s = e["scores"]
        values = [s["blosum_penalty"], s["hydro_norm"], s["charge_norm"], s["vol_norm"], s["grantham_norm"]]
        values_closed = values + [values[0]]
        fig.add_trace(go.Scatterpolar(
            r=values_closed, theta=categories_closed, fill="toself",
            name=e["label"], line_color=e["color"], opacity=0.65,
        ))
    fig.update_layout(
        polar=dict(
            bgcolor="rgba(0,0,0,0)",
            radialaxis=dict(visible=True, range=[0, 1], gridcolor=_BORDER_COLOR, linecolor=_BORDER_COLOR),
            angularaxis=dict(gridcolor=_BORDER_COLOR, linecolor=_BORDER_COLOR, tickfont=dict(size=11, family=_FONT)),
        ),
        showlegend=len(entries) > 1,
        legend=dict(font=dict(color=_TEXT_COLOR, size=11)),
    )
    return _base_layout(fig, height=380)


def make_waterfall_chart(scores: dict) -> go.Figure:
    """Waterfall chart depicting weight contribution to total MSI score."""
    metrics = ["BLOSUM62", "Hydrophobicity", "Charge Shift", "Volume Delta", "Grantham Dist"]
    contributions = [
        round(0.35 * scores["blosum_penalty"] * 100, 1),
        round(0.25 * scores["hydro_norm"] * 100, 1),
        round(0.20 * scores["charge_norm"] * 100, 1),
        round(0.10 * scores["vol_norm"] * 100, 1),
        round(0.10 * scores["grantham_norm"] * 100, 1),
    ]

    fig = go.Figure(go.Bar(
        x=metrics, y=contributions,
        marker_color=[COLOR_PRIMARY, "#3B82F6", "#8B5CF6", "#EC4899", "#F59E0B"],
        text=[f"+{v}" for v in contributions], textposition="auto",
        textfont=dict(family=_MONO, color="#000000", size=12),
    ))
    fig.update_layout(
        xaxis=dict(gridcolor=_BORDER_COLOR),
        yaxis=dict(title="MSI Point Contribution", gridcolor=_BORDER_COLOR, range=[0, 40]),
    )
    return _base_layout(fig, height=340)


def make_property_bar_chart(wt_aa: str, mut_aa: str) -> go.Figure:
    wt, mut = AA_PROPERTIES[wt_aa], AA_PROPERTIES[mut_aa]
    categories = ["Hydrophobicity", "Charge", "Volume (A^3 / 10)"]
    wt_vals = [wt["hydrophobicity"], wt["charge"], wt["volume"] / 10.0]
    mut_vals = [mut["hydrophobicity"], mut["charge"], mut["volume"] / 10.0]

    fig = go.Figure()
    fig.add_trace(go.Bar(name=f"Wildtype ({wt_aa})", x=categories, y=wt_vals, marker_color="#3B82F6"))
    fig.add_trace(go.Bar(name=f"Mutant ({mut_aa})", x=categories, y=mut_vals, marker_color=COLOR_PATHOGENIC))
    fig.update_layout(
        barmode="group",
        legend=dict(font=dict(color=_TEXT_COLOR, size=11)),
        yaxis=dict(gridcolor=_BORDER_COLOR, zerolinecolor="#374151"),
        xaxis=dict(gridcolor=_BORDER_COLOR),
    )
    return _base_layout(fig, height=340)


def make_distribution_histogram(all_scores: list, current_msi: float, color: str) -> go.Figure:
    fig = go.Figure()
    fig.add_trace(go.Histogram(
        x=all_scores, nbinsx=35, marker_color="#1E293B",
        marker_line_color="#334155", marker_line_width=1,
        name="Substitution Background",
    ))
    fig.add_vline(
        x=current_msi, line_width=2.5, line_dash="dash", line_color=color,
        annotation_text=f"Current Variant: {current_msi}",
        annotation_font_color=color, annotation_position="top left",
        annotation_font=dict(family=_MONO, size=12),
    )
    fig.update_layout(
        showlegend=False,
        xaxis=dict(title="Mutation Severity Index (MSI)", gridcolor=_BORDER_COLOR),
        yaxis=dict(title="Substitution Count", gridcolor=_BORDER_COLOR),
    )
    return _base_layout(fig, height=350)


def make_blosum_heatmap(wt_aa: str, mut_aa: str) -> go.Figure:
    fig = go.Figure(go.Heatmap(
        z=_BLOSUM62_ROWS, x=_AA_ORDER, y=_AA_ORDER,
        colorscale=[[0, "#7F1D1D"], [0.45, "#111827"], [1, "#06B6D4"]],
        colorbar=dict(title=dict(text="Score", font=dict(color=_TEXT_COLOR, size=11)), tickfont=dict(color=_TEXT_COLOR)),
        hovertemplate="WT: %{y} -> MUT: %{x} | Score: %{z}<extra></extra>",
    ))
    fig.add_trace(go.Scatter(
        x=[mut_aa], y=[wt_aa], mode="markers",
        marker=dict(symbol="square-open", size=22, line=dict(color="#FFFFFF", width=2.5)),
        showlegend=False, hoverinfo="skip",
    ))
    fig.update_layout(
        xaxis=dict(title="Mutant Residue", side="top", gridcolor=_BORDER_COLOR),
        yaxis=dict(title="Wildtype Residue", autorange="reversed", gridcolor=_BORDER_COLOR),
    )
    return _base_layout(fig, height=500)


def make_protein_schematic(gene: str, position: int, wt_aa: str, mut_aa: str) -> go.Figure:
    info = GENE_INFO[gene]
    length = info["protein_length"]
    domains = info["domains"]
    palette = ["#06B6D4", "#3B82F6", "#8B5CF6", "#EC4899", "#10B981", "#F59E0B"]

    fig = go.Figure()
    # Backbone track
    fig.add_shape(type="rect", x0=1, x1=length, y0=0.4, y1=0.6,
                  fillcolor="#1E293B", line=dict(color="#334155", width=1), layer="below")

    for i, d in enumerate(domains):
        fig.add_shape(type="rect", x0=d["start"], x1=d["end"], y0=0.25, y1=0.75,
                      fillcolor=palette[i % len(palette)], opacity=0.8,
                      line=dict(color="#FFFFFF", width=0.5))
        label_y = 0.95 if i % 2 == 0 else 1.15
        fig.add_annotation(x=(d["start"] + d["end"]) / 2, y=label_y, text=d["name"],
                            showarrow=False, font=dict(size=10, color=_TEXT_MUTED, family=_FONT))

    # Variant Marker
    fig.add_trace(go.Scatter(
        x=[position], y=[0.5], mode="markers+text",
        marker=dict(symbol="diamond", size=16, color=COLOR_PATHOGENIC, line=dict(color="#FFFFFF", width=1.5)),
        text=[f"p.{wt_aa}{position}{mut_aa}"], textposition="bottom center",
        textfont=dict(color=COLOR_PATHOGENIC, size=12, family=_MONO),
        showlegend=False,
    ))

    fig.update_xaxes(range=[-length * 0.03, length * 1.03], title=f"Residue Position (1 - {length} aa)",
                      gridcolor=_BORDER_COLOR, zeroline=False)
    fig.update_yaxes(range=[0, 1.35], visible=False)
    return _base_layout(fig, height=260)

# ------------------------------------------------------------------------------
# 4. CUSTOM STYLING (ENTERPRISE DARK THEME)
# ------------------------------------------------------------------------------

CUSTOM_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap');

/* Hide standard Streamlit header & decoration elements */
header[data-testid="stHeader"] { visibility: hidden; height: 0; }
#MainMenu { visibility: hidden; }
footer { visibility: hidden; }

/* Global Dark Theme Settings */
html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, sans-serif;
    background-color: #0B0F17 !important;
    color: #F3F4F6;
}

.stApp {
    background-color: #0B0F17;
}

/* Typography Overrides */
h1, h2, h3, h4, h5, h6 {
    font-family: 'Inter', sans-serif !important;
    font-weight: 700 !important;
    letter-spacing: -0.02em !important;
    color: #FFFFFF !important;
}

/* Sidebar Customization */
section[data-testid="stSidebar"] {
    background-color: #111827 !important;
    border-right: 1px solid #1F2937 !important;
}

/* Metric Display Overrides */
[data-testid="stMetricValue"] {
    font-family: 'JetBrains Mono', monospace !important;
    color: #06B6D4 !important;
    font-weight: 600 !important;
}
[data-testid="stMetricLabel"] {
    font-family: 'Inter', sans-serif !important;
    color: #9CA3AF !important;
    font-size: 0.8rem !important;
    text-transform: uppercase;
    letter-spacing: 0.05em;
}

/* Clean Professional Buttons */
.stButton > button {
    font-family: 'Inter', sans-serif;
    font-weight: 600;
    background-color: #06B6D4;
    color: #0B0F17;
    border: none;
    border-radius: 4px;
    padding: 0.5rem 1rem;
    transition: all 0.15s ease-in-out;
}
.stButton > button:hover {
    background-color: #22D3EE;
    color: #0B0F17;
    box-shadow: 0 0 12px rgba(6, 182, 212, 0.35);
}

/* Modern Minimalist Tabs */
button[data-baseweb="tab"] {
    font-family: 'Inter', sans-serif;
    font-weight: 600;
    font-size: 0.85rem;
    color: #9CA3AF;
    background-color: transparent;
    border: none !important;
    padding: 0.6rem 1rem;
}
button[data-baseweb="tab"][aria-selected="true"] {
    color: #06B6D4 !important;
    border-bottom: 2px solid #06B6D4 !important;
}

/* Custom Card Container */
.msi-card {
    background-color: #111827;
    border: 1px solid #1F2937;
    border-radius: 6px;
    padding: 1.0rem 1.2rem;
    margin-bottom: 0.8rem;
}
.msi-card-title {
    font-size: 0.75rem;
    color: #9CA3AF;
    text-transform: uppercase;
    letter-spacing: 0.08em;
    font-weight: 600;
    margin-bottom: 0.4rem;
}
.msi-card-value {
    font-size: 1.1rem;
    font-weight: 600;
    color: #F3F4F6;
}

/* System Status Pill */
.sys-badge {
    display: inline-block;
    padding: 0.2rem 0.6rem;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.72rem;
    border-radius: 3px;
    background: #1E293B;
    color: #06B6D4;
    border: 1px solid #334155;
    margin-bottom: 1rem;
}
</style>
"""

# ------------------------------------------------------------------------------
# 5. STREAMLIT APPLICATION CONTROLLER
# ------------------------------------------------------------------------------

st.set_page_config(
    page_title="Cancer Mutation Impact Predictor",
    layout="wide",
    initial_sidebar_state="expanded"
)
st.markdown(CUSTOM_CSS, unsafe_allow_html=True)

# Session State Initialization
if "history" not in st.session_state:
    st.session_state.history = []
if "current" not in st.session_state:
    st.session_state.current = None

# --- SIDEBAR CONTROLS ---
with st.sidebar:
    st.markdown("<div class='sys-badge'>GENOMIC ENGINE // v3.0</div>", unsafe_allow_html=True)
    st.markdown("<h3 style='margin-top:0;'>VARIANT CONFIGURATION</h3>", unsafe_allow_html=True)
    
    input_mode = st.radio("Input Method", ["Preset Driver Variants", "Custom Single Substitution"])

    preset_note = None
    if input_mode == "Preset Driver Variants":
        choice = st.selectbox("Select Benchmark Variant", list(PRESET_MUTATIONS.keys()))
        p = PRESET_MUTATIONS[choice]
        wt_aa, mut_aa, position, gene, preset_note = p["wt"], p["mut"], p["position"], p["gene"], p["note"]
    else:
        genes_list = list(GENE_INFO.keys()) + ["Custom / Unannotated"]
        gene_choice = st.selectbox("Target Gene", genes_list)
        gene = gene_choice if gene_choice != "Custom / Unannotated" else "Custom"
        
        aa_options = [f"{code} - {AA_PROPERTIES[code]['name']}" for code in AA_PROPERTIES]
        wt_sel = st.selectbox("Wildtype Amino Acid (WT)", aa_options, index=7) # Glycine default
        mut_sel = st.selectbox("Mutant Amino Acid (MUT)", aa_options, index=3) # Aspartate default
        position = st.number_input("Residue Position", min_value=1, value=12, step=1)
        wt_aa, mut_aa = wt_sel[0], mut_sel[0]

    st.markdown("---")
    execute_analysis = st.button("RUN PREDICTION ENGINE", use_container_width=True)
    
    st.markdown("---")
    st.markdown(f"<div style='font-size:0.8rem; color:#9CA3AF;'>Analyzed Variants in Session: <b>{len(st.session_state.history)}</b></div>", unsafe_allow_html=True)

# Process Trigger
if execute_analysis or st.session_state.current is None:
    if wt_aa == mut_aa:
        st.error("Synonymous Substitution Detected (WT == MUT). Mutation Severity Index is 0.0.")
    else:
        scores = compute_mutation_score(wt_aa, mut_aa)
        rec_id = len(st.session_state.history) + 1
        rec = {
            "id": rec_id, "gene": gene, "wt": wt_aa, "mut": mut_aa,
            "position": position, "scores": scores, "note": preset_note,
            "label": f"#{rec_id} {gene} {wt_aa}{position}{mut_aa}",
        }
        st.session_state.current = rec
        st.session_state.history.append(rec)

# --- MAIN DASHBOARD CONTENT ---
record = st.session_state.current

# Dashboard Header
st.markdown(
    """
    <div style="padding-bottom: 1.5rem; border-bottom: 1px solid #1F2937; margin-bottom: 1.5rem;">
        <h1 style="font-size: 2.0rem; margin-bottom: 0.2rem;">CANCER MUTATION IMPACT PREDICTOR</h1>
        <p style="color: #9CA3AF; font-size: 0.92rem; margin: 0;">
            Quantitative assessment of non-synonymous single amino acid substitutions combining BLOSUM62 conservation,
            Grantham evolutionary distance, and physicochemical property shifts.
        </p>
    </div>
    """,
    unsafe_allow_html=True,
)

if record:
    scores = record["scores"]
    wt_aa, mut_aa, position, gene = record["wt"], record["mut"], record["position"], record["gene"]
    verdict_text, verdict_color, verdict_key = get_verdict(scores["msi"])
    wt_name, mut_name = AA_PROPERTIES[wt_aa]["name"], AA_PROPERTIES[mut_aa]["name"]
    wt_cat, mut_cat, same_cat = classify_substitution(wt_aa, mut_aa)
    all_scores = get_all_substitution_scores()
    percentile = get_percentile(scores["msi"], all_scores)
    domain_text = get_domain_context(gene, position)

    # Variant Identifier Headline
    st.markdown(f"""
        <div style="display:flex; align-items:center; justify-content:space-between; background:#111827; padding:0.8rem 1.2rem; border-radius:6px; border:1px solid #1F2937; margin-bottom:1.2rem;">
            <div>
                <span style="font-family:'JetBrains Mono'; color:#06B6D4; font-weight:600; font-size:1.2rem;">{gene}</span>
                <span style="color:#F3F4F6; font-weight:600; font-size:1.2rem; margin-left:0.5rem;">p.{wt_name}{position}{mut_name} ({wt_aa}{position}{mut_aa})</span>
            </div>
            <div>
                <span style="background:{verdict_color}22; color:{verdict_color}; border:1px solid {verdict_color}66; font-family:'JetBrains Mono'; font-weight:600; padding:0.3rem 0.8rem; border-radius:4px; font-size:0.85rem;">
                    {verdict_text.upper()}
                </span>
            </div>
        </div>
    """, unsafe_allow_html=True)

    # Primary Navigation Tabs
    tabs = st.tabs([
        "OVERVIEW",
        "PHYSICOCHEMICAL PROFILE",
        "DISTRIBUTION & RANKING",
        "SUBSTITUTION MATRIX",
        "PROTEIN ARCHITECTURE",
        "SESSION COMPARISON",
        "METHODOLOGY"
    ])

    # --- TAB 1: OVERVIEW ---
    with tabs[0]:
        c1, c2 = st.columns([1.1, 1.9])
        with c1:
            st.plotly_chart(make_gauge_chart(scores["msi"], verdict_color), use_container_width=True)
            st.markdown(f"<div style='text-align:center; font-family:JetBrains Mono; color:{verdict_color}; font-size:1.1rem; margin-top:-1rem;'>MSI SCORE: {scores['msi']} / 100</div>", unsafe_allow_html=True)
        
        with c2:
            st.markdown("<div style='height:10px;'></div>", unsafe_allow_html=True)
            kpi1, kpi2, kpi3 = st.columns(3)
            with kpi1:
                st.markdown(f"<div class='msi-card'><div class='msi-card-title'>Domain Region</div><div class='msi-card-value'>{domain_text}</div></div>", unsafe_allow_html=True)
            with kpi2:
                st.markdown(f"<div class='msi-card'><div class='msi-card-title'>Class Transition</div><div class='msi-card-value'>{'Internal Class' if same_cat else 'Class Change'}</div></div>", unsafe_allow_html=True)
            with kpi3:
                st.markdown(f"<div class='msi-card'><div class='msi-card-title'>Percentile Rank</div><div class='msi-card-value'>Top {round(100-percentile,1)}%</div></div>", unsafe_allow_html=True)

            if gene in GENE_INFO:
                st.markdown(f"""
                    <div class='msi-card'>
                        <div class='msi-card-title'>Gene Annotations - {gene}</div>
                        <div style='font-size:0.88rem; color:#E5E7EB;'>
                            <b>{GENE_INFO[gene]['full_name']}</b><br>{GENE_INFO[gene]['function']}
                            {f'<br><br><span style="color:#06B6D4;">Clinical Note: ' + record['note'] + '</span>' if record['note'] else ''}
                        </div>
                    </div>
                """, unsafe_allow_html=True)

        st.markdown("---")
        exp1, exp2 = st.columns(2)
        with exp1:
            st.download_button(
                "EXPORT REPORT (TEXT)",
                data=generate_report_text(record),
                file_name=f"{gene}_{wt_aa}{position}{mut_aa}_report.txt",
                mime="text/plain",
                use_container_width=True
            )
        with exp2:
            st.download_button(
                "EXPORT REPORT (JSON)",
                data=generate_report_json(record),
                file_name=f"{gene}_{wt_aa}{position}{mut_aa}_report.json",
                mime="application/json",
                use_container_width=True
            )

    # --- TAB 2: PHYSICOCHEMICAL PROFILE ---
    with tabs[1]:
        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>COMPONENT WEIGHT CONTRIBUTIONS</div>", unsafe_allow_html=True)
            st.plotly_chart(make_waterfall_chart(scores), use_container_width=True)
        with col_b:
            st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>RAW PROPERTY SHIFTS (WT vs MUT)</div>", unsafe_allow_html=True)
            st.plotly_chart(make_property_bar_chart(wt_aa, mut_aa), use_container_width=True)

        st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>MULTI-PARAMETRIC SEVERITY RADAR</div>", unsafe_allow_html=True)
        st.plotly_chart(make_radar_chart([{"label": f"{wt_aa}->{mut_aa}", "scores": scores, "color": COLOR_PRIMARY}]), use_container_width=True)

    # --- TAB 3: DISTRIBUTION & RANKING ---
    with tabs[2]:
        st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>POSITION IN ALL 380 POSSIBLE SUBSTITUTIONS</div>", unsafe_allow_html=True)
        st.plotly_chart(make_distribution_histogram(all_scores, scores["msi"], verdict_color), use_container_width=True)
        st.markdown(f"""
            <div class='msi-card'>
                <div class='msi-card-title'>Statistical Severity Context</div>
                <div style='font-size:0.9rem; color:#E5E7EB;'>
                    The variant <b>p.{wt_aa}{position}{mut_aa}</b> yields an MSI score of <b>{scores['msi']}</b>, placing it higher in predicted severity than 
                    <b>{percentile}%</b> of all 380 non-synonymous single amino acid substitutions derived from the standard 20 amino acids.
                </div>
            </div>
        """, unsafe_allow_html=True)

    # --- TAB 4: SUBSTITUTION MATRIX ---
    with tabs[3]:
        st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>BLOSUM62 CONSERVATION MATRIX LOCATION</div>", unsafe_allow_html=True)
        st.plotly_chart(make_blosum_heatmap(wt_aa, mut_aa), use_container_width=True)

    # --- TAB 5: PROTEIN ARCHITECTURE ---
    with tabs[4]:
        if gene in GENE_INFO:
            st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>PROTEIN DOMAIN MAP & MUTATION LOCATION</div>", unsafe_allow_html=True)
            st.plotly_chart(make_protein_schematic(gene, position, wt_aa, mut_aa), use_container_width=True)
        else:
            st.info("Domain schematic mapping is available for annotated system genes (KRAS, TP53, EGFR, BRAF, PIK3CA, PTEN).")

    # --- TAB 6: SESSION COMPARISON ---
    with tabs[5]:
        if not st.session_state.history:
            st.info("No variants analyzed in the current session.")
        else:
            table_data = [
                {
                    "Session ID": h["id"],
                    "Gene": h["gene"],
                    "Mutation": f"{h['wt']}{h['position']}{h['mut']}",
                    "MSI Score": h["scores"]["msi"],
                    "Grantham Dist": h["scores"]["grantham_dist"],
                    "Verdict": get_verdict(h["scores"]["msi"])[0],
                }
                for h in st.session_state.history
            ]
            st.dataframe(table_data, use_container_width=True, hide_index=True)

            if len(st.session_state.history) >= 2:
                st.markdown("---")
                st.markdown("<div style='font-size:0.85rem; font-weight:600; color:#9CA3AF;'>OVERLAY RADAR COMPARISON</div>", unsafe_allow_html=True)
                labels = [h["label"] for h in st.session_state.history]
                selected = st.multiselect("Select Variants to Compare", labels, default=labels[-2:])
                if selected:
                    palette = ["#06B6D4", "#EF4444", "#F59E0B", "#10B981", "#8B5CF6"]
                    radar_entries = [
                        {
                            "label": h["label"],
                            "scores": h["scores"],
                            "color": palette[idx % len(palette)]
                        }
                        for idx, h in enumerate(st.session_state.history) if h["label"] in selected
                    ]
                    st.plotly_chart(make_radar_chart(radar_entries), use_container_width=True)

            if st.button("CLEAR SESSION HISTORY"):
                st.session_state.history = []
                st.session_state.current = None
                st.rerun()

    # --- TAB 7: METHODOLOGY ---
    with tabs[6]:
        st.markdown("""
        ### COMPUTATIONAL METHODOLOGY & ALGORITHMIC FORMULATION

        The **Mutation Severity Index (MSI)** is a deterministic heuristic scoring algorithm designed to quantify the disruption of single amino acid substitutions.

        #### Mathematical Model
        $$
        MSI = 100 \\times \\left[ 0.35 \\cdot P_{\\text{BLOSUM}} + 0.25 \\cdot H_{\\text{norm}} + 0.20 \\cdot C_{\\text{norm}} + 0.10 \\cdot V_{\\text{norm}} + 0.10 \\cdot G_{\\text{norm}} \\right]
        $$

        Where:
        * **$P_{\\text{BLOSUM}}$**: Normalized BLOSUM62 penalty derived as $\\max(0, (11 - S) / 15)$.
        * **$H_{\\text{norm}}$**: Normalized Kyte-Doolittle hydrophobicity shift absolute magnitude.
        * **$C_{\\text{norm}}$**: Formal net charge alteration ($|\\Delta C| / 2$).
        * **$V_{\\text{norm}}$**: Normalized van der Waals residue volume difference ($|\\Delta V| / 167.7$).
        * **$G_{\\text{norm}}$**: Grantham distance metric ($D_{\\text{Grantham}} / 215.0$).

        #### Primary Reference Data Sources
        1. **BLOSUM62 Matrix**: Henikoff S, Henikoff JG. *Amino acid substitution matrices from protein blocks.* PNAS, 1992.
        2. **Kyte-Doolittle Scale**: Kyte J, Doolittle RF. *A simple method for displaying the hydropathic character of a protein.* J Mol Biol, 1982.
        3. **Grantham Distance**: Grantham R. *Amino acid difference.* Science, 1974.
        4. **Residue Volumes**: Zamyatnin AA. *Protein volume in solution.* Prog Biophys Mol Biol, 1972.

        *Disclaimer: This software is intended for research and educational purposes only. It is not intended for clinical diagnosis or therapeutic decision-making.*
        """)

# --- FOOTER ---
st.markdown("---")
st.markdown(
    "<div style='text-align: center; color: #6B7280; font-size: 0.78rem; font-family: Inter;'>"
    "Cancer Mutation Impact Predictor &bull; Enterprise Bioinformatics Suite &bull; Non-Clinical Research Tool"
    "</div>",
    unsafe_allow_html=True
)