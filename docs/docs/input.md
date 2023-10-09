# Input data for PathIntegrate

## Omics data
Each multi-omics dataset should be in the form of a pandas DataFrame with samples as rows and molecules as columns. The index of the DataFrame should be the sample IDs, and the columns should be the molecule IDs. 

The molecule IDs should match those of the desired pathway database (i.e. ChEBI IDs, UniProt IDs, and ENSEMBL genes for Reactome; and KEGG IDs for KEGG). The values in the DataFrame should be the molecule abundances for each sample. 



| sample_id   | 1372     | 16610    | 72665    | 30915    | 37373    | Group     |
| ----------- | -------- | -------- | -------- | -------- | -------- | --------- |
| INCOV092-BL | 1.541009 | 1.228611 | 1.224076 | 1.962028 | 0.652984 | COVID     |
| INCOV107-BL | 0.910486 | 2.169111 | 2.819585 | 1.234384 | 1.453066 | COVID     |
| INCOV020-BL | 0.831297 | 0.23298       | 2.126393 | 0.861793 | 2.877589 | COVID     |
| INCOV035-BL | 1.862011 | 0.792962 | 1.434183 | 1.223473 | 0.706152 | COVID     |
| INCOV122-BL | 1.416927 | 2.493762 | 1.77004  | 0.888144 | 0.693444 | Non-covid |
| INCOV084-BL | 1.622171 | 1.021112 | 2.323956 | 0.573877 | 0.764003 | Non-covid |
| INCOV086-BL | 1.610941 | 1.205343 | 0.83498  | 2.600065 | 1.700068 | Non-covid |
| INCOV133-BL | 0.83727  | 2.144127 |  1.24222        | 1.035411 | 2.037335 | Non-covid |

## Pathways 
Pathways should be in the form of a pandas DataFrame with pathways as rows and molecules as columns. The index of the DataFrame should be the pathway IDs, and the columns should be the molecule IDs. The first column should be the pathway names or descriptions. 

Pathways can be from any pathway database, but the molecule IDs should match those of the omics data.

Each pathway can either contain molecules from a single omics, or a combination of omics.

|   | Pathway_name | 0      | 1      | 2      | 3      | 4      |
| - | ------------ | ------ | ------ | ------ | ------ | ------ |
| 1 | Pathway_52   | Q13554 | P61289 | P05114 | P62081 | P54760 |
| 2 | Pathway_53   | Q9Y243 | P17252 |        |        |        |
| 3 | Pathway_54   | 16708  | P06732 | P61289 | O00220 | O75914 |
| 4 | Pathway_55   | O15264 | P25786 |        |        |        |
| 5 | Pathway_56   | P07858 | P62979 | Q9Y625 | P14778 | P12314 |
| 6 | Pathway_57   | P18510 | P15260 | Q13557 | P32942 | P04818 |
| 7 | Pathway_58   | P00738 | P37023 | P01588 | P63098 | P05362 |
| 8 | Pathway_59   | P52798 | P15498 |        |        |        |