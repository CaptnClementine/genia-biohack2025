import streamlit as st
import pandas as pd
import sqlite3
from typing import Optional

st.set_page_config(page_title="GIA", layout="centered")
st.title("Gene Insights & Analysis")

st.markdown(
    """
    <style>
    [data-testid="stSidebar"] > div:first-child {
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
    }
    </style>
    """,
    unsafe_allow_html=True
)

@st.cache_resource(show_spinner=False)
def get_conn() -> Optional[sqlite3.Connection]:
    try:
        # Allow usage across Streamlit threads
        conn = sqlite3.connect("genes.db", check_same_thread=False)
        conn.row_factory = sqlite3.Row
        return conn
    except Exception as e:
        st.error(f"Could not connect to SQLite database: {e}")
        return None

@st.cache_data(show_spinner=False)
def load_gene_table() -> Optional[pd.DataFrame]:
    conn = get_conn()
    if not conn:
        return None
    try:
        df = pd.read_sql_query(
            "SELECT hgnc_symbol, ensembl_gene_id, description, chromosome_name FROM gene", conn
        )
        return df
    except Exception as e:
        st.error(f"Failed to load gene data from database: {e}")
        return None

@st.cache_data(show_spinner=False)
def load_consensus_for_ensembl(ensembl_ids: list[str]) -> Optional[pd.DataFrame]:
    if not ensembl_ids:
        return pd.DataFrame(columns=["Gene", "Tissue", "nTPM"])
    conn = get_conn()
    if not conn:
        return None
    try:
        placeholders = ",".join(["?"] * len(ensembl_ids))
        query = f"""
            SELECT g.ensembl_gene_id AS Gene, e.tissue AS Tissue, e.nTPM
            FROM protein_atlas_expression e
            JOIN gene g ON g.id = e.gene_id
            WHERE g.ensembl_gene_id IN ({placeholders})
        """
        df = pd.read_sql_query(query, conn, params=ensembl_ids)
        return df
    except Exception as e:
        st.error(f"Failed to load consensus data: {e}")
        return None

# --- New: CDoseMap loader ---
@st.cache_data(show_spinner=False)
def load_cdosemap_for_ensembl(ensembl_ids: list[str]) -> Optional[pd.DataFrame]:
    if not ensembl_ids:
        return pd.DataFrame(columns=["cnv_type","cytoband","size","discovery_sig","known_gd","gnomad_constrained_genes","db_source","article"])    
    conn = get_conn()
    if not conn:
        return None
    try:
        placeholders = ",".join(["?"] * len(ensembl_ids))
        q = f"""
            SELECT 
                   c.cnv_type,
                   c.cytoband,
                   c.size,
                   c.discovery_sig,
                   c.known_gd,
                   c.gnomad_constrained_genes,
                   c.db_source,
                   c.article
            FROM CDoseMap c
            JOIN gene g ON g.id = c.gene_id
            WHERE g.ensembl_gene_id IN ({placeholders})
        """
        return pd.read_sql_query(q, conn, params=ensembl_ids).dropna().drop_duplicates()
    except Exception as e:
        st.error(f"Failed to load CDoseMap data: {e}")
        return None


# Replace top-level UI with tabbed layout
gene_tab, about_tab = st.tabs(["Gene Search", "About"])

with gene_tab:
    df = load_gene_table()

    if df is not None and not df.empty:
        unique_vals = sorted(pd.unique(df['hgnc_symbol'].dropna().astype(str).str.strip())) if 'hgnc_symbol' in df.columns else []
        value_set_lower = {v.lower() for v in unique_vals}
        ensembl_vals = (
            sorted(pd.unique(df['ensembl_gene_id'].dropna().astype(str).str.strip()))
            if 'ensembl_gene_id' in df.columns else []
        )
        ensembl_set_lower = {v.lower() for v in ensembl_vals}

        # Move controls to sidebar
        with st.sidebar:
            st.subheader("Gene Search")
            manual = st.text_input(
                "Enter gene (HGNC symbol or Ensembl ID)",
                placeholder="e.g., BRCA1 or ENSG00000141510",
                key="manual_input",
            )
            selected = st.selectbox(
                "Select HGNC symbol",
                options=["-- Choose --"] + unique_vals,
                index=0,
                key="select_input",
            )

        # Decide which value to use (prefer manual input if provided)
        chosen = manual.strip() if manual and manual.strip() else (selected if selected != "-- Choose --" else "")

        # Main area: show results
        if chosen:
            chosen_l = chosen.casefold()
            in_hgnc = chosen_l in value_set_lower
            in_ensembl = chosen_l in ensembl_set_lower if ensembl_vals else False
            if in_hgnc or in_ensembl:
                st.header(chosen)
                mask = False
                if 'hgnc_symbol' in df.columns:
                    mask = df['hgnc_symbol'].astype(str).str.casefold().eq(chosen_l)
                if 'ensembl_gene_id' in df.columns:
                    mask = mask | df['ensembl_gene_id'].astype(str).str.casefold().eq(chosen_l)
                sel_rows = df[mask]
                display_cols = ['ensembl_gene_id', 'description', 'hgnc_symbol', 'chromosome_name']
                available_cols = [col for col in display_cols if col in sel_rows.columns]
                if available_cols:
                    sel_rows = sel_rows[available_cols]

                for index, row in sel_rows.iterrows():
                    for col in sel_rows.columns:
                        st.markdown(f"**{col}**: {row[col]}")
                        st.markdown("")

                # Consensus data from DB
                with st.expander("Consensus RNA data"):
                    if 'ensembl_gene_id' in sel_rows.columns:
                        selected_ensembl_ids = sel_rows['ensembl_gene_id'].dropna().astype(str).str.strip().tolist()
                        consensus_df = load_consensus_for_ensembl(selected_ensembl_ids)
                        if consensus_df is None or consensus_df.empty:
                            st.info("No consensus expression data found for the selected gene(s).")
                        else:
                            st.markdown("[Consensus RNA data from the Human Protein Atlas](https://www.proteinatlas.org/humanproteome/tissue/data#consensus_tissues_rna)")
                            st.markdown("The consensus normalized expression ('nTPM') value is calculated as the maximum nTPM value for each gene in the two data sources.")
                            consensus_filtered = consensus_df[['Tissue', 'nTPM']] if all(c in consensus_df.columns for c in ['Tissue','nTPM']) else consensus_df
                            st.dataframe(consensus_filtered, use_container_width=True, hide_index=True)

                # CDoseMap data
                with st.expander("CDoseMap data"):
                    if 'ensembl_gene_id' in sel_rows.columns:
                        cd_ids = sel_rows['ensembl_gene_id'].dropna().astype(str).str.strip().tolist()
                        cd_df = load_cdosemap_for_ensembl(cd_ids)
                        if len(cd_df) < 1:
                            st.markdown("No CDoseMap data found for the selected gene")
                            st.image("tg_image_3445653099.jpeg", width = 200)
                        else:
                            show_cols = [c for c in [
                                'cnv_type','cytoband','size','discovery_sig','known_gd',
                                'gnomad_constrained_genes','db_source','article'
                            ] if c in cd_df.columns]
                            st.markdown("Data from [the paper](https://www.cell.com/cell/fulltext/S0092-8674(22)00788-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867422007887%3Fshowall%3Dtrue#mmc1)")
                            st.markdown("Data were provided with gene symbols only. Do not extrapolate to specific Ensembl IDs")
                            st.dataframe(cd_df[show_cols], use_container_width=True, hide_index=True)
            else:
                st.warning("We don't have such data")
                suggestions = [v for v in unique_vals if v.lower().startswith(chosen.lower())]
                if ensembl_vals:
                    suggestions += [v for v in ensembl_vals if v.lower().startswith(chosen.lower())]
                if suggestions:
                    st.info("Did you mean: " + ", ".join(suggestions[:5]))

with about_tab:
    try:
        with open("APP_README.md", "r", encoding="utf-8") as f:
            content = f.read()
        st.markdown(content)
    except FileNotFoundError:
        st.warning("APP_README.md not found in the project root.")
    except Exception as e:
        st.error(f"Could not read APP_README.md: {e}")
