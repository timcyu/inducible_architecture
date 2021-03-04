## Description of processed datasets

#### `fLP3.loop_MOPS_glu_expression.txt`
Expression data for the 615 lacUV5 variants described in Figure 1. The `name` column can be parsed to obtain the TF identity, offset, and proximal site identity. The `RNA_exp_12_IPTG` and `RNA_exp_12` columns correspond to induced and uninduced expression respectively.

#### `LacZ_variant_exp.txt`
Expression data for the 8269 lacUV5 variants (Pcombo, Pmultiple, Pcore, Psteric) described in Figures 2-5. The `name` column can parsed to obtain the architecture label and identity of each promoter part. The `normalized_RNA_exp_Induced_12`, `normalized_RNA_exp_UnInduced_12`, and `ratio` columns correspond to induced expression, uninduced expression, and fold-change respectively.

Please note: Promoters with 'DISTAL' in their `name` are actually the Pcore promoters described. Additionally, promoters with 'HYBRID' in their `name` are not described at all.

#### `modelExp.txt`
Predicted Pcombo expression data by thermodynamic model for reproducing Figure 3D. 
