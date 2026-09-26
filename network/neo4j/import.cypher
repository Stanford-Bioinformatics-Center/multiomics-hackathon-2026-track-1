// =====================================================================================================
// network/neo4j/import.cypher — LOAD THE EXERCISE NETWORKS INTO NEO4J (5.x)
// =====================================================================================================
//
// PURPOSE: builds the graph described in network/neo4j/README.md from the CSV files written by
//   export_neo4j.R. Safe to re-run: every node and relationship is MERGEd on its key, so re-importing an
//   updated export updates properties instead of duplicating anything.
// HOW TO RUN: copy the files from $HACK_OUT/neo4j_import/ into Neo4j's import folder, then either paste this
//   file into Neo4j Browser (enable multi-statement queries) or run
//     cypher-shell -u neo4j -p <password> -f network/neo4j/import.cypher
// CONVENTIONS: empty CSV cells become absent properties; ";"-separated text becomes number lists; weights
//   are w_EE / w_RE / w_diff = w_EE - w_RE (positive = higher in endurance). Relationships between two nodes
//   of the same kind are stored once (direction has no meaning; query them undirected: -[r]-).
// =====================================================================================================

// ---- keys (one node per key; also makes lookups fast) ----------------------------------------------------
CREATE CONSTRAINT gene_key IF NOT EXISTS FOR (g:Gene) REQUIRE g.entrez IS UNIQUE;
CREATE CONSTRAINT metabolite_key IF NOT EXISTS FOR (m:Metabolite) REQUIRE m.name IS UNIQUE;
CREATE CONSTRAINT class_key IF NOT EXISTS FOR (c:MetaboliteClass) REQUIRE c.name IS UNIQUE;
CREATE CONSTRAINT contrast_key IF NOT EXISTS FOR (k:Contrast) REQUIRE k.id IS UNIQUE;
CREATE CONSTRAINT normalisation_key IF NOT EXISTS FOR (n:Normalisation) REQUIRE n.ome IS UNIQUE;
CREATE CONSTRAINT dataset_key IF NOT EXISTS FOR (d:Dataset) REQUIRE d.name IS UNIQUE;
CREATE INDEX gene_symbol IF NOT EXISTS FOR (g:Gene) ON (g.symbol);

// ---- gene / protein nodes -----------------------------------------------------------------------------------
LOAD CSV WITH HEADERS FROM 'file:///nodes_gene.csv' AS r
MERGE (g:Gene {entrez: r.entrez_gene})
SET g.symbol = r.symbol, g.uniprot = r.uniprot,
    g.in_string = (r.in_string = 'TRUE'), g.string_degree = toInteger(r.string_degree),
    g.string_component = toInteger(r.string_component), g.string_component_size = toInteger(r.string_component_size),
    g.mean_response_EE = toFloat(r.mean_response_EE), g.mean_response_RE = toFloat(r.mean_response_RE),
    g.joint_strength_EE = toFloat(r.joint_strength_EE), g.joint_strength_RE = toFloat(r.joint_strength_RE),
    g.joint_degree = toInteger(r.joint_degree),
    g.x_gene = toFloat(r.x_gene), g.y_gene = toFloat(r.y_gene),
    g.x_joint = toFloat(r.x_joint), g.y_joint = toFloat(r.y_joint),
    g.x_joint14 = toFloat(r.x_joint14), g.y_joint14 = toFloat(r.y_joint14),
    g.feature_id_adipose_rna = r.feature_id_adipose_rna, g.feature_id_adipose_prot = r.feature_id_adipose_prot,
    g.feature_id_blood_rna = r.feature_id_blood_rna, g.feature_id_blood_prot = r.feature_id_blood_prot,
    g.feature_id_muscle_rna = r.feature_id_muscle_rna, g.feature_id_muscle_prot = r.feature_id_muscle_prot,
    g.embedding_EE = [x IN split(r.embedding_EE, ';') | toFloat(x)],
    g.embedding_RE = [x IN split(r.embedding_RE, ';') | toFloat(x)],
    g.embedding_dims = split(r.embedding_dims, ';');

// ---- metabolite nodes -------------------------------------------------------------------------------------
LOAD CSV WITH HEADERS FROM 'file:///nodes_metabolite.csv' AS r
MERGE (m:Metabolite {name: r.name})
SET m.refmet_name = r.refmet_name, m.refmet_id = r.refmet_id, m.super_class = r.super_class, m.main_class = r.main_class,
    m.pubchem_cid = r.pubchem_cid, m.inchi_key = r.inchi_key, m.chebi_id = r.chebi_id, m.chebi_all = r.chebi_all,
    m.lookup_status = r.lookup_status,
    m.platform_adipose = r.platform_adipose, m.platform_blood = r.platform_blood, m.platform_muscle = r.platform_muscle,
    m.mean_response_EE = toFloat(r.mean_response_EE), m.mean_response_RE = toFloat(r.mean_response_RE),
    m.joint_strength_EE = toFloat(r.joint_strength_EE), m.joint_strength_RE = toFloat(r.joint_strength_RE),
    m.joint_degree = toInteger(r.joint_degree),
    m.x_metab = toFloat(r.x_metab), m.y_metab = toFloat(r.y_metab),
    m.x_joint = toFloat(r.x_joint), m.y_joint = toFloat(r.y_joint),
    m.x_joint14 = toFloat(r.x_joint14), m.y_joint14 = toFloat(r.y_joint14),
    m.embedding_EE = [x IN split(r.embedding_EE, ';') | toFloat(x)],
    m.embedding_RE = [x IN split(r.embedding_RE, ';') | toFloat(x)],
    m.embedding_dims = split(r.embedding_dims, ';'),
    m.embedding_doubled_EE = [x IN split(r.embedding_doubled_EE, ';') | toFloat(x)],
    m.embedding_doubled_RE = [x IN split(r.embedding_doubled_RE, ';') | toFloat(x)],
    m.embedding_doubled_dims = split(r.embedding_doubled_dims, ';');

// ---- metabolite classes, contrasts, normalisation, dataset ------------------------------------------------
LOAD CSV WITH HEADERS FROM 'file:///nodes_class.csv' AS r
MERGE (c:MetaboliteClass {name: r.name})
SET c.level = 'RefMet super class', c.n_metabolites = toInteger(r.n_metabolites), c.n_in_joint_network = toInteger(r.n_in_joint_network);

LOAD CSV WITH HEADERS FROM 'file:///nodes_contrast.csv' AS r
MERGE (k:Contrast {id: r.contrast_id})
SET k.arm = r.arm, k.arm_label = r.arm_label, k.tissue = r.tissue, k.ome = r.ome, k.time = r.time, k.time_h = toFloat(r.time_h);

LOAD CSV WITH HEADERS FROM 'file:///nodes_normalisation.csv' AS r
MERGE (n:Normalisation {ome: r.ome})
SET n.divisor_max_abs_logFC = toFloat(r.divisor_max_abs_logFC), n.set_by = r.set_by, n.n_values = toInteger(r.n_values);

LOAD CSV WITH HEADERS FROM 'file:///nodes_dataset.csv' AS r
MERGE (d:Dataset {name: r.name})
SET d.exported_at = r.exported_at, d.code_commit = r.code_commit, d.source_package = r.source_package,
    d.string_cutoff = toInteger(r.string_cutoff), d.rhea_release = r.rhea_release, d.normalisation = r.normalisation;

// ---- network edges ------------------------------------------------------------------------------------------
// gene - gene: STRING combined score >= 700 (step 2), weights per arm (step 3)
LOAD CSV WITH HEADERS FROM 'file:///rel_string.csv' AS r
MATCH (a:Gene {entrez: r.entrez_a}), (b:Gene {entrez: r.entrez_b})
MERGE (a)-[e:STRING_INTERACTS]->(b)
SET e.combined_score = toInteger(r.combined_score),
    e.w_EE = toFloat(r.w_EE), e.w_RE = toFloat(r.w_RE), e.w_diff = toFloat(r.w_diff),
    e.sig_EE = toFloat(r.sig_EE), e.sig_RE = toFloat(r.sig_RE), e.cos_EE = toFloat(r.cos_EE), e.cos_RE = toFloat(r.cos_RE);

// metabolite - metabolite: shared or STRING-linked Rhea enzyme + same super class (step 6)
LOAD CSV WITH HEADERS FROM 'file:///rel_class_link.csv' AS r
MATCH (a:Metabolite {name: r.metabolite_a}), (b:Metabolite {name: r.metabolite_b})
MERGE (a)-[e:CLASS_LINK]->(b)
SET e.class = r.class, e.class_level = r.class_level, e.link_type = r.link_type,
    e.n_shared_proteins = toInteger(r.n_shared_proteins),
    e.shared_proteins = CASE WHEN r.shared_proteins IS NULL THEN [] ELSE split(r.shared_proteins, ';') END,
    e.string_protein_pairs = CASE WHEN r.string_protein_pairs IS NULL THEN [] ELSE split(r.string_protein_pairs, ';') END,
    e.w_EE = toFloat(r.w_EE), e.w_RE = toFloat(r.w_RE), e.w_diff = toFloat(r.w_diff),
    e.sig_EE = toFloat(r.sig_EE), e.sig_RE = toFloat(r.sig_RE);

// metabolite - protein: Rhea (step 5), weights = doubled-metabolite-embedding dot products (step 14)
LOAD CSV WITH HEADERS FROM 'file:///rel_rhea.csv' AS r
MATCH (m:Metabolite {name: r.metabolite}), (g:Gene {entrez: r.entrez_gene})
MERGE (m)-[e:RHEA_ENZYME]->(g)
SET e.n_reactions = toInteger(r.n_reactions), e.example_reactions = split(r.example_reactions, ';'),
    e.matched_via = r.matched_via, e.w_EE = toFloat(r.w_EE), e.w_RE = toFloat(r.w_RE), e.w_diff = toFloat(r.w_diff);

// metabolite -> class
LOAD CSV WITH HEADERS FROM 'file:///rel_in_class.csv' AS r
MATCH (m:Metabolite {name: r.metabolite}), (c:MetaboliteClass {name: r.super_class})
MERGE (m)-[:IN_CLASS]->(c);

// ---- exercise responses (one relationship per node x contrast) ----------------------------------------------
LOAD CSV WITH HEADERS FROM 'file:///rel_gene_response.csv' AS r
CALL {
  WITH r
  MATCH (g:Gene {entrez: r.entrez_gene}), (k:Contrast {id: r.contrast_id})
  MERGE (g)-[x:RESPONDS_IN]->(k)
  SET x.logFC = toFloat(r.logFC), x.normalised = toFloat(r.normalised), x.se = toFloat(r.se)
} IN TRANSACTIONS OF 5000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///rel_metabolite_response.csv' AS r
CALL {
  WITH r
  MATCH (m:Metabolite {name: r.metabolite}), (k:Contrast {id: r.contrast_id})
  MERGE (m)-[x:RESPONDS_IN]->(k)
  SET x.logFC = toFloat(r.logFC), x.normalised = toFloat(r.normalised), x.se = toFloat(r.se)
} IN TRANSACTIONS OF 5000 ROWS;

// contrasts -> the normalisation their values were divided by
MATCH (k:Contrast), (n:Normalisation) WHERE n.ome = k.ome
MERGE (k)-[:NORMALISED_BY]->(n);

// ---- check: counts per label and relationship type (compare with README "Expected counts") ----------------
MATCH (n) RETURN labels(n)[0] AS label, count(*) AS n ORDER BY label;
MATCH ()-[r]->() RETURN type(r) AS relationship, count(*) AS n ORDER BY relationship;
