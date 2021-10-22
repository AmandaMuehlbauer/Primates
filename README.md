# Primates
These are the scripts and processed data associated with the publication “Interspecies variation in hominid gut microbiota controls host gene regulation.”

The file subread_DESeq2-data.RDATA is the DESeq2 object. 

The taxonomic results from the HUMAnN2 pipeline are Merged_Class.txt, Merged_Family.txt, Merged_Genus.txt, Merged_Order.txt, Merged_Phylum.txt and Merged_Species.txt. The pathway data is in PrimatePathAbundance.tsv. 

The general analyses are in FinalAnalyses_v09_github.R. This includes the framework of likelihood ratio tests used to classify each gene according to which hominid microbiomes it responds to. It also includes the basic DESeq2 models used throughout the paper. 

The rest of the scripts are for making the figures that appear in the paper. 
