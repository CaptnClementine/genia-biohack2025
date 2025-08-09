import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import pandas as pd

def create_violin_plot(gene_symbols_for_plot: list[str], plot_df: pd.DataFrame, stat: str, features=None):
    """
    Create a violin plot comparing expression levels for selected genes vs other genes.
    There should be 2 violins plot on the same plot
    1 - only data for gene_symbols_for_plot
    2 - data for all genes without data for gene_symbols_for_plot
    y - values of plot_df['stat']
    """
    #import pdb; pdb.set_trace()
    gene_col = "ensembl_gene_id"
    # Create two datasets: selected genes vs other genes
    selected_data = plot_df[plot_df[gene_col].isin(gene_symbols_for_plot)].copy()
    other_data = plot_df[~plot_df[gene_col].isin(gene_symbols_for_plot)].copy()
    
    # Add group labels
    selected_data['Group'] = 'Selected Genes'
    other_data['Group'] = 'Other Genes'
    
    # Combine data
    combined_data = pd.concat([selected_data, other_data], ignore_index=True)
    
    if combined_data.empty:
        st.info("No data available for plotting.")
        return
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Create violin plot with two groups
    sns.violinplot(
        data=combined_data,
        x="Group",
        y=stat,
        inner=None,
        cut=0,
        linewidth=0,
        palette=["lightblue", "lightcoral"],
        ax=ax
    )
    
    # Make violins transparent
    for coll in ax.collections:
        try:
            coll.set_alpha(0.35)
        except Exception:
            pass
    
    # Overlay boxplot
    sns.boxplot(
        data=combined_data,
        x="Group",
        y=stat,
        showcaps=True,
        width=0.3,
        whiskerprops={'linewidth': 1.1},
        medianprops={'color': 'black', 'linewidth': 1.2},
        showfliers=False,
        ax=ax,
        palette='Accent'
    )
    
  
    ax.set_xlabel("")
    ax.set_ylabel(stat)
    ax.set_title("Expression Comparison: Selected vs Other Genes")
    
    sns.despine(ax=ax)
    fig.tight_layout()
    st.pyplot(fig)