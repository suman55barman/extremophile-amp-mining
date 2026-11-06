import streamlit as st
import pandas as pd
import numpy as np
from pathlib import Path
import py3Dmol
import streamlit.components.v1 as components
import plotly.express as px
import plotly.graph_objects as go

# Page setup
st.set_page_config(
    page_title="Extremophile AMP Discovery Engine",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

BASE_DIR = Path(__file__).resolve().parent
LEADS_DIR = BASE_DIR / "data" / "05_final_leads"
REPORT_TSV = LEADS_DIR / "prioritized_leads_report.tsv"
PDBS_DIR = LEADS_DIR / "pdbs"

st.title("Extreophile AMP Mining Dashboard")
st.caption("Autonomous Bioprospecting, Safety Screening & Structural Validation")

if not REPORT_TSV.exists():
    st.error(f"Report file not found at {REPORT_TSV}. Run `master_amp_mining.py` first.")
    st.stop()

# Load data
df = pd.read_csv(REPORT_TSV, sep="\t")

# Top metrics Row
c1, c2, c3, c4 = st.columns(4)
c1.metric("Total Prioritized Leads", len(df))
c2.metric("Top Priority Score", f"{df['priority_score'].max():.1f}" if 'priority_score' in df.columns else "N/A")
c3.metric("Avg Hydrophobic Moment", f"{df['uH'].mean():.2f}" if 'uH' in df.columns else "N/A")
c4.metric("Avg Net Charge", f"{df['net_charge'].mean():.2f}" if 'net_charge' in df.columns else "N/A")

st.divider()

# Layout: Table on left, 3D structure viewer on right
col_table, col_viewer = st.columns([1.2, 1])

with col_table:
    st.subheader("Prioritized Candidate Leaderboard")
    st.dataframe(
        df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "sequence": st.column_config.TextColumn("Sequence", width="medium"),
            "priority_score": st.column_config.ProgressColumn("Score", min_value=0, max_value=100, format="%.1f")
        }
    )

with col_viewer:
    st.subheader("ESMFold 3D Structure Viewer")

    # List Available PDBs
    pdb_files = sorted(list(PDBS_DIR.glob("*.pdb")))
    if not pdb_files:
        st.info("No PDB files found in data/05_final_leads/pdbs/")
    else:
        selected_pdb_name = st.selectbox(
            "Select Candidate PDB to Inspect",
            options=[p.name for p in pdb_files]
        )
        selected_pdb_path = PDBS_DIR / selected_pdb_name

        # Color scheme selector
        color_scheme = st.radio("Color By", ["spectrum", "hydrophobicity", "chain"], horizontal=True)

        with open(selected_pdb_path, "r") as f:
            pdb_data = f.read()

        # Render 3D mol view
        view = py3Dmol.view(width=500, height=450)
        view.addModel(pdb_data, "pdb")
        view.setStyle({"cartoon": {"color": color_scheme}})
        view.zoomTo()
        components.html(view._make_html(), height=450, width=500)

        st.caption(f"Rendering: `{selected_pdb_name}` via py3Dmol")