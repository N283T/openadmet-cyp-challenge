SELECT act.activity_id, act.molregno, m.chembl_id AS molecule_chembl_id,
       s.canonical_smiles AS smiles, s.standard_inchi_key,
       t.chembl_id AS target_chembl_id, t.pref_name AS target_name, t.target_type,
       t.organism AS target_organism, a.assay_id, a.chembl_id AS assay_chembl_id,
       a.description AS assay_description, a.assay_type, a.confidence_score,
       a.assay_organism, a.src_id AS assay_src_id, a.src_assay_id,
       act.standard_type, act.standard_relation, act.standard_value::double precision,
       act.standard_units, act.pchembl_value::double precision,
       act.activity_comment, act.data_validity_comment, act.potential_duplicate,
       act.doc_id, d.chembl_id AS document_chembl_id, d.pubmed_id, d.doi, d.title,
       d.year, d.src_id AS document_src_id
FROM activities act
JOIN molecule_dictionary m ON m.molregno = act.molregno
LEFT JOIN compound_structures s ON s.molregno = act.molregno
JOIN assays a ON a.assay_id = act.assay_id
JOIN target_dictionary t ON t.tid = a.tid
LEFT JOIN docs d ON d.doc_id = act.doc_id
WHERE t.organism = 'Homo sapiens'
  AND t.pref_name IN ('Cytochrome P450 1A2', 'Cytochrome P450 2C9',
                      'Cytochrome P450 2D6', 'Cytochrome P450 3A4', 'Cytochrome P450 2C19')
  AND act.standard_type IN ('IC50', 'AC50', 'Potency', 'Ki')
ORDER BY act.activity_id
