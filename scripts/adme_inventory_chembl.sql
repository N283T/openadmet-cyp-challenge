SELECT a.chembl_id AS assay_chembl_id, a.description AS assay_description,
       a.assay_type, a.assay_organism, a.confidence_score, a.src_id,
       act.standard_type, act.standard_units, act.standard_relation,
       count(*) AS observation_rows, count(DISTINCT act.molregno) AS molecules,
       count(*) FILTER (WHERE act.data_validity_comment IS NULL
                       AND coalesce(act.potential_duplicate,0)=0) AS unflagged_rows,
       min(act.standard_value)::double precision AS minimum,
       max(act.standard_value)::double precision AS maximum
FROM assays a JOIN activities act ON act.assay_id=a.assay_id
WHERE a.src_id=27
GROUP BY a.chembl_id,a.description,a.assay_type,a.assay_organism,a.confidence_score,
         a.src_id,act.standard_type,act.standard_units,act.standard_relation
ORDER BY a.chembl_id,act.standard_type,act.standard_units,act.standard_relation;
