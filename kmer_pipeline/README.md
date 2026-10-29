### Pipeline for identifying clonal isolates using k-mer profiles ###

- v3scrapeSRA uses a .txt file containing a unique SRA or ENA identifier, after which the following structure produces the following files:
- subsampler.py: creates subsampled .fastq files of single or paired end data
- kmercounting.sh: produces k-mer files using Jellyfish per sample
- jellyfish_to_parquet.py: takes the jellyfish files produced by kmercounting.sh and converts them to .parquet files
- pairwisecomparison.py: uses the .parquet files produced by jellyfish_to_parquet.py to create Jaccard, weighted Jaccard and Dice-sorensen similarity matrices

All scripts can be automated using the sample_per_time_pipeline.py script
