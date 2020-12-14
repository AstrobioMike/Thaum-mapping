The main output files are in the 'combined-results' directory, which holds tables of combined results across all samples:

Gene-level stuff

	- combined-results/Combined-gene-coverages-and-KO-annots.tsv
		- coverage values for all genes across all samples (columns are: gene ID; KO ID; KO function name; then the samples)

	- combined-results/Combined-CPM-normalized-gene-coverages-and-KO-annots.tsv
		- coverage values for all genes across all samples normalized to coverage per million (CPM) – (each coverage value from the above table is divided by the sum of its column (sample), then multiplied by 1000000), so more comparable across samples, but keep in mind this is still normalized to total coverage, so the ultimate denominator here is dependent on the Thaum population that was detected (I think this often makes sense when we want to look at changes within/among the Thaum populations, but there is also a table with read recruitment percentages i'll mention below)


Genome-level stuff

	- combined-results/Combined-genome-coverages.tsv
		- average genome coverages across all samples

	- combined-results/Combined-CPM-normalized-genome-coverages.tsv
		-same deal as mentioned above, genome-level coverage here is normalized within-sample to be like a percent, but times 1 million instead of times 100 (same note as above too to keep in mind, that these CPM normalizations are of total coverage, and therefore are a window focused on the Thaum populations we detected in each sample)

	- combined-results/Combined-genome-detections.tsv
		- genome detections across all samples (fraction of total basepairs that recruited at least 1 read)

	- combined-results/Combined-genome-detection-filtered-coverages.tsv
		- average genome coverages across samples after filtering by detection (if less than 50% of a reference genome's basepairs recruited any reads, we set the coverage value for that genome in that sample to 0 in this file)

	- combined-results/Combined-CPM-normalized-genome-detection-filtered-coverages.tsv
		- same deal as discussed above, normalized within sample to coverage per million, but this took the detection-filtered table and normalized that


Sample recruitment overview

	- combined-results/Per-sample-mapping-overview.tsv
		- each row is a sample, columns include: number of fragments (counts pairs or single-end reads the same way, e.g. the "fragment" is the starting DNA molecule that would have become a single-end read or a paired-end read); if the sample was single-end or paired-end; and ther percentage of overall reads aligned as reported by bowtie2 for that sample (this is to all of our reference genomes, so if wanted, this can be how to get a pseudo-relative abundance of Thaum reads recovered from each sample)


Genes sequences, gff files, and KO annotations are all in 'genes-and-annotations'. 


Reference genome seqs and download scripts/logs are in 'ref-genomes'.


Workflow and code are in 'metagenome-dl-and-mapping-workflow'.


"reference-genomes.xlsx" holds the reference genome spreadsheet we worked with to identify our ref genome targets
"MG-samples.xlsx" holds the metagenomes spreadsheet we worked with to identify our sample targets
