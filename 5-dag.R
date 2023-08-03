# 6 DAGITY: 
library(dagitty)

dag <- dagitty('dag {
  bb="0,0,1,1"
  "bf/formula_ratio" [pos="0.637,0.567"]
  BM_collection_time [pos="0.616,0.772"]
  HMO_level [exposure,pos="0.197,0.855"]
  breastfeeding_frequency [pos="0.630,0.662"]
  diet [pos="0.061,0.749"]
  ethnicity [pos="0.109,0.659"]
  genetics [pos="0.061,0.855"]
  lactation_phase [pos="0.298,0.643"]
  miRNA_level [outcome,pos="0.818,0.855"]
  season [pos="0.193,0.642"]
  "bf/formula_ratio" -> HMO_level
  "bf/formula_ratio" -> miRNA_level
  BM_collection_time -> HMO_level
  BM_collection_time -> miRNA_level
  HMO_level -> miRNA_level
  breastfeeding_frequency -> HMO_level
  breastfeeding_frequency -> miRNA_level
  diet -> HMO_level
  ethnicity -> HMO_level
  genetics -> HMO_level
  lactation_phase -> HMO_level
  season -> HMO_level
}')
plot(dag)

dag2 <- dagitty('dag {
bb="0,0,1,1"
"bf/formula_ratio" [pos="0.637,0.567"]
BM_collection_time [pos="0.616,0.772"]
HMO_level [exposure,pos="0.197,0.855"]
breastfeeding_frequency [pos="0.630,0.662"]
diet [pos="0.061,0.749"]
ethnicity [pos="0.109,0.659"]
genetics [pos="0.061,0.855"]
lactation_phase [pos="0.298,0.643"]
miRNA_level [outcome,pos="0.818,0.855"]
rRNA_proportion [pos="0.748,0.514"]
season [pos="0.193,0.642"]
supernatant_volume [pos="0.847,0.511"]
"bf/formula_ratio" -> HMO_level
"bf/formula_ratio" -> miRNA_level
BM_collection_time -> HMO_level
BM_collection_time -> miRNA_level
HMO_level -> miRNA_level
breastfeeding_frequency -> HMO_level
breastfeeding_frequency -> miRNA_level
diet -> HMO_level
ethnicity -> HMO_level
genetics -> HMO_level
lactation_phase -> HMO_level
rRNA_proportion -> miRNA_level
season -> HMO_level
supernatant_volume -> miRNA_level
}')

