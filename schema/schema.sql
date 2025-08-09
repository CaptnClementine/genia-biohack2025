-- Gene database schema
-- This schema defines the structure for storing gene information and protein atlas expression data
DROP table if exists gene;
CREATE TABLE gene (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    hgnc_symbol TEXT ,
    ensembl_gene_id TEXT not null,
    chromosome_name TEXT,
    description TEXT,
    UNIQUE(ensembl_gene_id)
);
DROP table if exists protein_atlas_expression;

CREATE TABLE protein_atlas_expression (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    tissue TEXT NOT NULL,
    nTPM REAL,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE,
    UNIQUE(gene_id, tissue)
);

-- CDoseMap CNV data
DROP TABLE IF EXISTS CDoseMap;
CREATE TABLE CDoseMap (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gene_id INTEGER NOT NULL,
    cnv_type TEXT,
    cytoband TEXT,
    size INTEGER,
    discovery_sig TEXT,
    known_gd TEXT,
    gnomad_constrained_genes TEXT,
    db_source TEXT,
    article TEXT,
    FOREIGN KEY(gene_id) REFERENCES gene(id) ON DELETE CASCADE
);

-- Indexes for better query performance
CREATE INDEX idx_gene_hgnc_symbol ON gene(hgnc_symbol);
CREATE INDEX idx_gene_ensembl_gene_id ON gene(ensembl_gene_id);
CREATE INDEX idx_expr_gene ON protein_atlas_expression(gene_id);
CREATE INDEX idx_expr_tissue ON protein_atlas_expression(tissue);
CREATE INDEX idx_cdose_gene ON CDoseMap(gene_id);
CREATE INDEX idx_cdose_cnv_type ON CDoseMap(cnv_type);
CREATE INDEX idx_cdose_cytoband ON CDoseMap(cytoband);
