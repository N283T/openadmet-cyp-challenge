SELECT act.activity_id,act.molregno,m.chembl_id AS molecule_chembl_id,
 s.canonical_smiles AS raw_smiles,s.standard_inchi_key,
 a.chembl_id AS source_assay_id,a.description AS assay_description,a.assay_organism,
 a.assay_type,a.confidence_score,a.src_id,
 act.standard_type,act.standard_units,act.standard_relation,
 act.standard_value::double precision,act.data_validity_comment,act.potential_duplicate,
 act.activity_comment,act.doc_id,d.chembl_id AS document_chembl_id,d.doi,d.title
FROM assays a JOIN activities act ON act.assay_id=a.assay_id
JOIN molecule_dictionary m ON m.molregno=act.molregno
LEFT JOIN compound_structures s ON s.molregno=act.molregno
LEFT JOIN docs d ON d.doc_id=act.doc_id
WHERE a.src_id=27 AND a.chembl_id IN
 ('CHEMBL3301363','CHEMBL3301364','CHEMBL3301365','CHEMBL3301370','CHEMBL3301372')
ORDER BY act.activity_id;
