SELECT t.chembl_id AS target_chembl_id, t.pref_name, t.organism, t.target_type,
       a.assay_type, a.confidence_score, a.assay_organism,
       act.standard_type, act.standard_relation, act.standard_units,
       act.data_validity_comment, act.potential_duplicate,
       count(*) AS rows, count(DISTINCT act.molregno) AS molecules,
       count(DISTINCT act.assay_id) AS assays,
       count(act.pchembl_value) AS with_pchembl,
       count(act.standard_value) AS with_standard_value,
       min(act.pchembl_value) AS min_pchembl, max(act.pchembl_value) AS max_pchembl
FROM activities act
JOIN assays a ON a.assay_id = act.assay_id
JOIN target_dictionary t ON t.tid = a.tid
WHERE t.organism = 'Homo sapiens'
  AND t.pref_name IN ('Cytochrome P450 1A2', 'Cytochrome P450 2C9',
                      'Cytochrome P450 2D6', 'Cytochrome P450 3A4')
GROUP BY t.chembl_id, t.pref_name, t.organism, t.target_type,
         a.assay_type, a.confidence_score, a.assay_organism,
         act.standard_type, act.standard_relation, act.standard_units,
         act.data_validity_comment, act.potential_duplicate
ORDER BY t.pref_name, act.standard_type, act.standard_relation,
         a.assay_type, a.confidence_score
