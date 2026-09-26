// =====================================================================================================
// network/neo4j/queries.cypher — EXAMPLE QUERIES FOR THE VISUALISER (run after import.cypher)
// =====================================================================================================
// Each query is independent; run them one at a time in Neo4j Browser, or all with cypher-shell -f.
// Weights: w_EE, w_RE (per arm) and w_diff = w_EE - w_RE (positive = higher in endurance). Not tested
// statistically: treat differences as descriptive.
// =====================================================================================================

// 1. A gene and its neighbours in the joint network (STRING partners and Rhea metabolites), largest |w_diff| first.
MATCH (g:Gene {symbol: 'SRC'})-[e:STRING_INTERACTS|RHEA_ENZYME]-(n)
RETURN g.symbol AS gene, type(e) AS edge, coalesce(n.symbol, n.name) AS neighbour,
       round(e.w_EE, 4) AS w_EE, round(e.w_RE, 4) AS w_RE, round(e.w_diff, 4) AS w_diff
ORDER BY abs(e.w_diff) DESC;

// 2. The 10 edges of any type whose weight differs most between the arms.
MATCH (a)-[e:STRING_INTERACTS|CLASS_LINK|RHEA_ENZYME]-(b) WHERE elementId(a) < elementId(b)
RETURN coalesce(a.symbol, a.name) AS a, coalesce(b.symbol, b.name) AS b, type(e) AS edge, round(e.w_diff, 4) AS w_diff,
       CASE WHEN e.w_diff > 0 THEN 'higher in endurance' ELSE 'higher in resistance' END AS direction
ORDER BY abs(e.w_diff) DESC LIMIT 10;

// 3. One metabolite class with its enzymes: the class subgraph to draw (e.g. fatty acyls).
MATCH (c:MetaboliteClass {name: 'Fatty Acyls'})<-[:IN_CLASS]-(m:Metabolite)
OPTIONAL MATCH (m)-[r:RHEA_ENZYME]->(g:Gene)
RETURN m.name AS metabolite, collect(g.symbol) AS enzymes, round(m.mean_response_EE, 4) AS mean_EE, round(m.mean_response_RE, 4) AS mean_RE
ORDER BY size(enzymes) DESC;

// 4. A gene's exercise responses across tissues, omes and times, both arms (the data behind its embedding).
MATCH (g:Gene {symbol: 'SRC'})-[x:RESPONDS_IN]->(k:Contrast)
RETURN k.arm AS arm, k.tissue AS tissue, k.ome AS ome, k.time AS time, round(x.logFC, 3) AS logFC, round(x.normalised, 4) AS normalised, round(x.se, 3) AS se
ORDER BY arm, tissue, ome, k.time_h;

// 5. Check: a metabolite - protein weight is one dot product of the doubled metabolite embedding and the gene embedding.
MATCH (m:Metabolite)-[e:RHEA_ENZYME]->(g:Gene)
WITH m, e, g, reduce(s = 0.0, i IN range(0, size(g.embedding_EE) - 1) | s + m.embedding_doubled_EE[i] * g.embedding_EE[i]) AS dot_EE
RETURN count(*) AS rhea_edges, max(abs(dot_EE - e.w_EE)) AS largest_mismatch;

// 6. Check: a gene - gene weight is the dot product of the two genes' embeddings.
MATCH (a:Gene)-[e:STRING_INTERACTS]->(b:Gene)
WITH e, reduce(s = 0.0, i IN range(0, size(a.embedding_RE) - 1) | s + a.embedding_RE[i] * b.embedding_RE[i]) AS dot_RE
RETURN count(*) AS string_edges, max(abs(dot_RE - e.w_RE)) AS largest_mismatch;

// 7. Everything a visualiser needs for the joint network in one pull: nodes with layout + colours inputs, edges with weights.
MATCH (n) WHERE (n:Gene OR n:Metabolite) AND n.x_joint IS NOT NULL
WITH collect({id: coalesce(n.symbol, n.name), kind: labels(n)[0], class: n.super_class, x: n.x_joint, y: n.y_joint,
              mean_EE: n.mean_response_EE, mean_RE: n.mean_response_RE, strength_EE: n.joint_strength_EE, strength_RE: n.joint_strength_RE}) AS nodes
MATCH (a)-[e:STRING_INTERACTS|CLASS_LINK|RHEA_ENZYME]->(b)
RETURN size(nodes) AS n_nodes, count(e) AS n_edges, nodes[0] AS example_node,
       collect({from: coalesce(a.symbol, a.name), to: coalesce(b.symbol, b.name), type: type(e), w_EE: e.w_EE, w_RE: e.w_RE, w_diff: e.w_diff})[0] AS example_edge;

// 8. Provenance: which export and normalisation the database holds.
MATCH (d:Dataset) OPTIONAL MATCH (n:Normalisation)
RETURN d.exported_at AS exported_at, d.code_commit AS code_commit, collect(n.ome + ' / ' + toString(round(n.divisor_max_abs_logFC, 3)) + ' (' + n.set_by + ')') AS divisors;
