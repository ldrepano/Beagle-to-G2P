# Beagle-to-G2P

Code to facilitate the analysis of Base Editing data with the Genomics to Proteins (G2P) Portal. In collaboration with Ganna Reint (Genetic Perturbation Platform R&D, Broad Institute of MIT and Harvard). 

## generate-G2P-input-file

The python script in this folder associates guide-level base editing data with individual residues along a specified transcript for integration with the G2P portal. Guides with multiple residues in their editing window are associated with the residue in which they make the most severe mutation. For residues targeted by numerous guides, the residue is associated with the guides with the most extreme z-score. 

Usage:
```
python3 Beagle_to_G2P.py sample_Beagle_to_G2P_input.txt
```
Input file should specify:
1. Path to file containing Z-scores for each sgRNA tiling the gene. 
	- Can be excel, csv, or tab-delimited 
	- Must include following columns: sgRNA Target Sequence, Z-score
2. Path to Per-codon base edits Beagle output file 
3. A-G or C-T depending on whether the screen employed an A-based or C-based editor 
4. Select transcript to use
	- must be present in the "Target Transcript ID" column of codon-level Beagle output
	- ideally canonical MANE Select Transcript
5. Identifier to prepend to output file name 

Example input is provided in folder. 

Note that, in output, residues for which the strongest mutation possible is indicated but are lacking information on the "strongest Z-score sgRNA" are uniquely targeted by guides that introduce a more severe edit at neighboring codons. 

## visualize-G2P-results

The python script in this folder uses the **output** file from G2P (following the pipeline from generate-G2P-input-file) the generate an interactive plot. 

Usage:
```
python3 BE_G2P_output_visualization.py G2P_PPM1D_O15297_protein_features_CBE.csv
```
G2P_PPM1D_O15297_protein_features_CBE.csv represents a sample G2P output file and is included in the folder. 

Point shape is determined by the mutation type introduced by the sgRNA associated with the reported Z-score. Point color is determined by features introduced by any additional columns in the file. For generalizability, the plot projects all features that are present in the G2P output file. To select/deselect certain features for visualization, select them in the legend. 

**NOTE:** Edits are matched with ClinVar-annotated mutations on amino acid *position*, not identity. 
