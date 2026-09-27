"""Transparent protein-exerkine prioritization from pinned human MoTrPAC summaries.

No participant-level model is fit. Source BH is computed across complete assay
families by extract_exerkine_screen.R before gene matching; Python verifies plasma.
Run from the notebook, or directly to regenerate derived tables and provenance.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ['EE-CON', 'RE-CON']
FAMILY = ['tissue', 'assay', 'platform', 'contrast']
INPUTS = ['HUMAN_FEATURE_TO_GENE', 'BLOOD_PROT_OL_DA', 'BLOOD_TRNSCRPT_DA',
          'MUSCLE_TRNSCRPT_DA', 'ADIPOSE_TRNSCRPT_DA',
          'MUSCLE_PROT_PR_DA', 'ADIPOSE_PROT_PR_DA']
FOCAL = ['CD300LG', 'ANGPT2', 'WARS1', 'CCN1', 'CX3CL1',
         'FLT3LG', 'ANGPTL7', 'HGF', 'STC2', 'TIMP3']
SOURCES = {
    'motrpac': 'https://github.com/MoTrPAC/MotrpacHumanPreSuspensionAnalysis/tree/535b4044e7417413de471104c619120337602b77',
    'landscape': 'https://doi.org/10.64898/2026.02.27.702183',
    'CD300LG': 'https://doi.org/10.7554/eLife.96535.3',
    'ANGPT2': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC5391203/',
    'WARS1': 'https://doi.org/10.1016/j.celrep.2022.111905',
    'FLT3LG': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9999360/',
    'ANGPTL7': 'https://doi.org/10.1371/journal.pone.0173024',
    'atlas': 'https://exerkineatlas.org/',
}


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    assert len(p) and np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p, kind='stable')
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty(len(p))
    out[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1)
    return out


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_inputs(root=ROOT):
    source = root / 'motrpac-exploration'
    prior = json.loads((source / 'table3/input_provenance.json').read_text())
    hashes = {name + '.rda': sha256(source / (name + '.rda')) for name in INPUTS}
    for name, value in hashes.items():
        assert value == prior['input_sha256'][name], f'Input changed: {name}'
    return hashes


def annotate_time(df):
    """Post-bout minutes; negative values only order during-EE samples.

    During-20 and during-40 correspond to -20/0 relative to the end of
    the 40-minute endurance bout. No during-resistance samples exist.
    """
    out = df.copy()
    minute = {'during_20_min': -20, 'during_40_min': 0, 'post_10_min': 10,
              'post_24_hr': 1440}
    out['minutes_from_bout_end'] = out.Timepoint.map(minute)
    early = out.Timepoint.eq('post_15_30_45_min')
    late = out.Timepoint.eq('post_3.5_4_hr')
    out.loc[early, 'minutes_from_bout_end'] = out.loc[early, 'tissue'].map(
        {'blood': 30, 'muscle': 15, 'adipose': 45})
    out.loc[late, 'minutes_from_bout_end'] = out.loc[late, 'tissue'].map(
        {'blood': 210, 'muscle': 210, 'adipose': 240})
    assert out.minutes_from_bout_end.notna().all()
    def label(row):
        if row.Timepoint.startswith('during'):
            return row.Timepoint.replace('_', ' ')
        m = row.minutes_from_bout_end
        return f'post {int(m)} min' if m < 60 else f'post {m / 60:g} h'
    out['sample_time'] = out.apply(label, axis=1)
    return out


def read_annotations(reference):
    """Use structured locations, not free-text keyword guesses of secretion.

    A signal peptide also occurs on membrane proteins. Its presence alone is
    insufficient. Missing secretion annotation is not evidence of non-secretion.
    Nonclassical secretion literature is assessed separately in curated notes.
    """
    rows = []
    for filename in sorted(reference.glob('uniprot_[0-9]*.json')):
        for entry in json.loads(filename.read_text())['results']:
            locations = []
            for comment in entry.get('comments', []):
                if comment['commentType'] == 'SUBCELLULAR LOCATION':
                    locations += [x['location']['value']
                                  for x in comment.get('subcellularLocations', [])]
            rows.append({
                'uniprot': entry['primaryAccession'],
                'uniprot_gene': ';'.join(g.get('geneName', {}).get('value', '')
                                       for g in entry.get('genes', [])),
                'entry_type': entry['entryType'],
                'annotation_date': entry.get('entryAudit', {}).get('lastAnnotationUpdateDate'),
                'locations': '; '.join(sorted(set(locations))),
                'location_has_secreted': any(x.startswith('Secreted') for x in locations),
                'location_has_membrane': any('membrane' in x.lower() for x in locations),
                'signal_peptide': any(f['type'] == 'Signal' for f in entry.get('features', [])),
                'uniprot_url': 'https://www.uniprot.org/uniprotkb/' + entry['primaryAccession'] + '/entry',
            })
    ann = pd.DataFrame(rows)
    assert not ann.uniprot.duplicated().any()
    requests = json.loads((reference / 'uniprot_requests.json').read_text())
    requested = {a for request in requests for a in request['accessions']}
    assert requested == set(ann.uniprot), 'Incomplete UniProt download'
    return ann


def join_symbols(series):
    return ';'.join(sorted(set(series.dropna().astype(str))))


def read_published_candidates(reference):
    path = reference / 'table_s8.xlsx'
    assert hashlib.md5(path.read_bytes()).hexdigest() == 'beaf62ed40e4618b6655078fef49d310'
    table = pd.read_excel(path, sheet_name='Exerkine_Candidate_List')
    # Explicit aliases relevant to our focal genes; WAS is a different gene.
    table['canonical_gene'] = table.gene_symbol.replace(
        {'WARS': 'WARS1', 'CYR61': 'CCN1', 'WISP1': 'CCN4'})
    return table


def screen(plasma, tissue, alpha):
    """Heuristic intersection, NOT a gene-level hypothesis test or FDR estimate."""
    p = plasma[(plasma.adj_p_value < alpha) & (plasma.logFC > 0)]
    t = tissue[(tissue.adj_p_value < alpha) & (tissue.logFC > 0)]
    cols = ['feature_id', 'Timepoint', 'sample_time', 'minutes_from_bout_end',
            'logFC', 'p_value', 'adj_p_value', 'contrast', 'tissue', 'assay']
    pairs = p[['gene_symbol', 'contrast_category'] + cols].merge(
        t[['gene_symbol', 'contrast_category'] + cols],
        on=['gene_symbol', 'contrast_category'], suffixes=('_plasma', '_tissue'))
    pairs['tissue_at_or_before_plasma'] = (
        pairs.minutes_from_bout_end_tissue <= pairs.minutes_from_bout_end_plasma)
    summary = pairs.groupby('gene_symbol').agg(
        min_plasma_q=('adj_p_value_plasma', 'min'),
        min_tissue_q=('adj_p_value_tissue', 'min'),
        matched_modes=('contrast_category', join_symbols),
        supporting_tissues=('tissue_tissue', join_symbols),
        supporting_assays=('assay_tissue', join_symbols),
        any_tissue_at_or_before_plasma=('tissue_at_or_before_plasma', 'any'),
        significant_time_pairs=('contrast_category', 'size')).reset_index()
    return p, pairs, summary.sort_values('min_plasma_q')


def curation():
    """Explicit manual prioritization; no numeric score or inferred causality."""
    rows = [
        ('CD300LG', 'Priority: diabetes follow-up', 'Membrane protein; circulating molecular form unresolved',
         'Reported as an exercise-responsive candidate in a 2024 human training study; associations with insulin sensitivity and genetic analyses support glucose-homeostasis follow-up, not proof of therapeutic benefit.',
         'Identify the measured circulating fragment/full-length protein by orthogonal assays; then test effects on endothelial function and glucose handling.', SOURCES['CD300LG']),
        ('ANGPT2', 'Priority: vascular remodeling', 'Secreted ligand in UniProt',
         'Already highlighted in MoTrPAC Figure 8. Adipose-specific mouse work links ANGPT2-driven angiogenesis to metabolism; acute human plasma increases do not establish benefit.',
         'Test whether post-exercise plasma changes endothelial TIE2 signaling and barrier function, and whether ANGPT2 depletion/add-back changes that response.', SOURCES['landscape'] + ' | ' + SOURCES['ANGPT2']),
        ('WARS1', 'Priority: exploratory immune signaling', 'Nonclassical secretion demonstrated outside exercise',
         'Cell experiments demonstrate soluble and vesicular WARS1 release and innate inflammatory activity. Absent from the checked MoTrPAC Table S8, including alias WARS. Exercise-related release route, active species, and disease effect remain unestablished; broader novelty is not claimed.',
         'Resolve soluble versus vesicular WARS1 and cell-injury markers; test WARS1 depletion/add-back in an innate-immune reporter assay.', SOURCES['WARS1']),
        ('CCN1', 'Previously analyzed benchmark', 'Secreted matrix-associated protein',
         'Positive control from the existing CCN1 notebook, also highlighted in MoTrPAC Figure 8.',
         'Use to benchmark the screen; do not count as a new candidate.', SOURCES['landscape']),
        ('CX3CL1', 'Previously analyzed benchmark', 'Membrane and soluble forms',
         'Positive control from the existing fractalkine notebook; listed in cached Exerkine Atlas.',
         'Use to benchmark the screen; distinguish membrane-bound protein from soluble fractalkine.', SOURCES['atlas']),
        ('FLT3LG', 'Secondary: plasma response only', 'Secreted and membrane forms annotated',
         'Strong plasma response, also described in a prior human exercise study. No positive muscle/adipose RNA or total-protein hit under the primary rule.',
         'Investigate immune-cell sources and independent circulating-protein validation.', SOURCES['FLT3LG']),
        ('ANGPTL7', 'Secondary: plasma response only', 'Secreted extracellular protein annotated',
         'Acute endurance plasma increase here; a prior obesity study reported decreases after a training program. Different populations and acute versus chronic sampling prevent direct equivalence.',
         'Test whether an acute pulse and longer-term baseline change differ in paired samples.', SOURCES['ANGPTL7']),
        ('HGF', 'Exploratory: tissue misses 0.05', 'Established secreted growth factor; structured-location annotation alone is incomplete',
         'Plasma passes 0.05; corresponding positive muscle RNA q is just above 0.05. Listed only in the threshold-sensitivity results.',
         'Replicate the tissue response before promoting to the primary cross-tissue list.', 'https://www.uniprot.org/uniprotkb/P14210/entry'),
        ('STC2', 'Exploratory: plasma misses 0.05', 'Secreted protein annotated',
         'Strong muscle RNA increase but plasma q is above 0.05. Not a primary plasma hit.',
         'Use an independent plasma assay or cohort; do not lower the cutoff after observing the result.', 'https://www.uniprot.org/uniprotkb/O76061/entry'),
        ('TIMP3', 'Exploratory: plasma misses 0.05', 'Secreted extracellular-matrix protein annotated',
         'Already highlighted in MoTrPAC Figure 8 using a broader screen; this pinned analysis has plasma q above 0.05.',
         'Check independent circulating-protein evidence; retain the preset 0.05 primary threshold.', SOURCES['landscape']),
    ]
    return pd.DataFrame(rows, columns=['gene_symbol', 'priority', 'release_interpretation',
                                       'evidence_and_limit', 'testable_follow_up', 'source_urls'])


def build_screen(root=ROOT):
    root = Path(root)
    out = root / 'results/exerkine_screen'
    ref = out / 'reference_data'
    hashes = verify_inputs(root)
    audit = pd.read_csv(out / 'bh_audit.csv')
    assert (audit.max_BH_difference < 1e-10).all()
    data = annotate_time(pd.read_csv(out / 'mapped_source_results.csv.gz', low_memory=False))
    full = annotate_time(pd.read_csv(out / 'plasma_full_background.csv'))
    full['bh_python'] = full.groupby(FAMILY, observed=True, dropna=False).p_value.transform(bh_adjust)
    np.testing.assert_allclose(full.bh_python, full.adj_p_value, rtol=1e-10, atol=1e-12)
    assert not full.duplicated(FAMILY + ['feature_id']).any()
    assert full.groupby(FAMILY, observed=True).size().eq(1417).all()
    valid = full.gene_symbol.notna() & full.ambiguous_gene.eq(False)
    full['pooled_primary_plasma_bh'] = np.nan
    primary_mask = full.contrast_category.isin(PRIMARY)
    # A sensitivity analysis across all 1,417 x 10 primary plasma tests.
    # This is a different family, not the original source correction.
    full.loc[primary_mask, 'pooled_primary_plasma_bh'] = bh_adjust(full.loc[primary_mask, 'p_value'])
    plasma = full[valid & primary_mask].copy()
    tissue = data[data.tissue.isin(['muscle', 'adipose']) &
                  data.assay.isin(['transcript-rna-seq', 'prot-pr']) &
                  data.contrast_category.isin(PRIMARY)].copy()
    ann = read_annotations(ref)
    ann.to_csv(out / 'uniprot_annotations.csv', index=False)
    # Aggregate annotations over the measured UniProt accessions for each gene.
    gene_ann = []
    for gene, group in plasma.groupby('gene_symbol'):
        ids = sorted({a for val in group.uniprot.dropna() for a in val.split(';')})
        a = ann[ann.uniprot.isin(ids)]
        gene_ann.append({'gene_symbol': gene, 'uniprot': ';'.join(ids),
                         'annotation_available': len(a) == len(ids) and len(ids) > 0,
                         'locations': join_symbols(a.locations),
                         'location_has_secreted': bool(a.location_has_secreted.any()),
                         'location_has_membrane': bool(a.location_has_membrane.any()),
                         'signal_peptide': bool(a.signal_peptide.any())})
    gene_ann = pd.DataFrame(gene_ann)
    atlas = pd.read_csv(root / 'motrpac-exploration/atlas_protein_gene_mapping.csv')
    gene_ann['in_cached_atlas'] = gene_ann.gene_symbol.isin(atlas.gene_symbol)
    published = read_published_candidates(ref)
    gene_ann['in_published_table_s8'] = gene_ann.gene_symbol.isin(published.canonical_gene)
    score = published.groupby('canonical_gene').Extracellular_score.max()
    gene_ann['published_COMPARTMENTS_score'] = gene_ann.gene_symbol.map(score)
    published.to_csv(out / 'published_table_s8_extracted.csv', index=False)
    notes = curation()
    p05, pairs05, summary05 = screen(plasma, tissue, .05)
    p10, pairs10, summary10 = screen(plasma, tissue, .10)
    assert gene_ann.set_index('gene_symbol').loc[p10.gene_symbol.unique(), 'annotation_available'].all()
    summaries = {}
    for name, p, pairs, summary in [('primary', p05, pairs05, summary05),
                                     ('exploratory_010', p10, pairs10, summary10)]:
        s = summary.merge(gene_ann, on='gene_symbol', validate='one_to_one').merge(
            notes, on='gene_symbol', how='left', validate='one_to_one')
        s['priority'] = s.priority.fillna('Hold: extracellular signaling evidence requires review')
        s.to_csv(out / f'{name}_same_mode_candidates.csv', index=False)
        pairs.to_csv(out / f'{name}_evidence_pairs.csv', index=False)
        p.to_csv(out / f'{name}_plasma_positive_hits.csv', index=False)
        summaries[name] = s
    assert {'CCN1', 'CX3CL1', 'CD300LG', 'ANGPT2', 'WARS1'} <= set(summary05.gene_symbol)
    # All plasma-up genes stay visible, including those without tissue overlap.
    all_up = p05.groupby('gene_symbol').agg(
        min_plasma_q=('adj_p_value', 'min'), plasma_modes=('contrast_category', join_symbols),
        n_positive_plasma_tests=('feature_id', 'size')).reset_index().merge(gene_ann, on='gene_symbol')
    all_up['passes_same_mode_tissue_rule'] = all_up.gene_symbol.isin(summary05.gene_symbol)
    all_up = all_up.sort_values('min_plasma_q')
    all_up.to_csv(out / 'all_primary_plasma_up_genes.csv', index=False)
    notes.to_csv(out / 'candidate_curation.csv', index=False)
    full.to_csv(out / 'plasma_background_validated.csv', index=False)
    # Keep nonsignificant estimates and direct modality tests for all focal genes.
    focal = data[data.gene_symbol.isin(FOCAL)].copy()
    focal = focal.merge(full[['feature_id', 'contrast', 'pooled_primary_plasma_bh']],
                        on=['feature_id', 'contrast'], how='left', validate='many_to_one')
    focal.to_csv(out / 'focal_genes_all_source_contrasts.csv', index=False)
    full[full.gene_symbol.isin(FOCAL) & primary_mask].to_csv(out / 'focal_plasma_all_primary.csv', index=False)
    # Ensure the previous single-gene results are reproduced exactly.
    benchmark_counts = {}
    for gene, path in [('CCN1', 'results/ccn1/ccn1_plasma_primary.csv'),
                       ('CX3CL1', 'results/fractalkine/cx3cl1_plasma_differential.csv')]:
        old = pd.read_csv(root / path)
        check = plasma[plasma.gene_symbol.eq(gene)].merge(old, on=['feature_id', 'contrast'], suffixes=('_new', '_old'), validate='one_to_one')
        assert len(check) == len(old) == 10
        for col in ['logFC', 'p_value', 'adj_p_value']:
            np.testing.assert_allclose(check[col + '_new'], check[col + '_old'], rtol=1e-10, atol=1e-12)
        benchmark_counts[gene] = len(check)
    counts = {
        'source_BH_families_validated': int(len(audit)),
        'source_feature_contrasts_validated': int(audit.n_features.sum()),
        'max_BH_difference': float(audit.max_BH_difference.max()),
        'plasma_assay_features': int(full.feature_id.nunique()),
        'plasma_primary_tests': int(primary_mask.sum()),
        'plasma_features_unmapped_or_ambiguous': int(full[~valid].feature_id.nunique()),
        'primary_plasma_up_genes': int(p05.gene_symbol.nunique()),
        'primary_same_mode_genes': int(len(summary05)),
        'primary_explicit_secreted_location_genes': int(summaries['primary'].location_has_secreted.sum()),
        'primary_any_tissue_at_or_before_plasma': int(summary05.any_tissue_at_or_before_plasma.sum()),
        'exploratory_010_plasma_up_genes': int(p10.gene_symbol.nunique()),
        'exploratory_010_same_mode_genes': int(len(summary10)),
        'pooled_BH_positive_plasma_genes': int(full[primary_mask & valid & full.logFC.gt(0) & full.pooled_primary_plasma_bh.lt(.05)].gene_symbol.nunique()),
        'uniprot_accessions_annotated': int(len(ann)),
        'published_table_s8_genes': int(published.canonical_gene.nunique()),
        'primary_genes_in_published_table_s8': int(summaries['primary'].in_published_table_s8.sum()),
        'primary_genes_absent_published_table_s8': sorted(set(summary05.gene_symbol) - set(published.canonical_gene)),
        'benchmarks_reproduced': benchmark_counts,
    }
    (out / 'screen_counts.json').write_text(json.dumps(counts, indent=2) + '\n')
    reference_hashes = {str(p.relative_to(root)): sha256(p) for p in sorted(ref.iterdir())
                        if p.is_file() and p.suffix in {'.json', '.xml', '.xlsx'}}
    provenance = {
        'analysis_date': '2026-09-26', 'dataset': 'human-precovid-sed-adu c2.0',
        'package_version': '2.0.8', 'motrpac_commit': '535b4044e7417413de471104c619120337602b77',
        'input_sha256': hashes, 'reference_sha256': reference_hashes,
        'source_urls': SOURCES,
        'published_candidates_source': 'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13184684/supplementaryFiles',
        'published_candidates_file': 'media-8.xlsx / Exerkine_Candidate_List; original MD5 verified against article XML',
        'raw_p_values': 'Supplied by MoTrPAC; no participant-level model refit.',
        'BH_family': FAMILY,
        'primary_rule': 'Plasma logFC>0 and source BH<0.05; positive muscle/adipose RNA or total protein source BH<0.05 for same gene and exercise mode, at any sampled time.',
        'threshold_sensitivity': 'Separate exploratory screen at 0.10; not primary discoveries.',
        'pooled_plasma_sensitivity': 'BH over all 14170 primary plasma feature-comparison raw p values; different family from source BH.',
        'annotation_rule': 'UniProt structured location starts with Secreted; membrane/signal peptide flags separate; published nonclassical/membrane release evidence curated separately.',
        'scope': 'Protein exerkines measured by Olink; not metabolites, unmeasured proteins, or tissue-restricted paracrine signals.',
        'counts': counts,
    }
    (out / 'input_provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    return {'full': full, 'data': data, 'plasma': plasma, 'tissue': tissue,
            'primary': summaries['primary'], 'exploratory': summaries['exploratory_010'],
            'pairs': pairs05, 'all_up': all_up, 'curation': notes,
            'focal': focal, 'audit': audit, 'counts': counts}


if __name__ == '__main__':
    result = build_screen()
    print(json.dumps(result['counts'], indent=2))
    print(result['primary'][['gene_symbol', 'min_plasma_q', 'min_tissue_q', 'matched_modes', 'priority']].to_string(index=False))
