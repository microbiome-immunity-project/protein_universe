Below we are providing workflows to search the MIP database against a protein sequence or structure or function of interest.

# SEQUENCE SEARCH

FASTA files (see `metadata/fasta` folder) contain all raw/curated MIP sequences for Rosetta and DMPFold models and may be directly used for sequence search/clustering etc. Sequence header ID corresponds to the MIP ID whereas header description to the NCBI sequence identifier, e.g.:

`>MIP_00004878 NC_013093.1|896|8-61`

`MIP_00004878` # MIP ID
`NC_013093.1` # NCBI ID
`896` # genome ID
`8-61` # sequence range

In order to facilitate sequence search we generated BLAST databases using `ncbi-blast-2.13.0+`:

```
makeblastdb -in rosetta_curated.fasta -dbtype prot -parse_seqids -out mip_rosetta_curated_db
```

and similarly for other databases. Sequence search may be performed e.g. using `blastp`:

```
blastp -query test.fasta -db mip_rosetta_curated_db
```

Note 1: The databases (sequences, structures and functions) are separated by quality metrics such that Rosetta and DMPFold datasets are in different datasets, which is why `rosetta_curated.fasta` contains the sequences and MIP IDs that belong to the Rosetta models in the curated dataset. The files for `rosetta_raw.fasta`, `dmpfold_curated.fasta`, and `dmpfold_faw.fasta` are named accordingly.

Note 2: Taxonomic information is not explicitly included in the databases but using descriptions in the FASTA headers it is easy to identify taxonomy.

# STRUCTURAL COMPARISON

There are a variety of ways one might use to find the structural relationships between a given protein and the MIP models. The method described below uses [TMAlign++ from the Zhang group](https://zhanggroup.org/TM-align/).

## Unpacking the models

`TMAlign++` requires the PDB files to be untar'd and uncompressed. It cannot operate on the PDB files in the MIP dataset without first uncompressing them. It is recommended to do this in parallel. This can be accomplished easily using standard command line tools like `find` and `xargs`.

Command lines to uncompress the Rosetta folded, high-quality model PDB files using 16 processes are below. The Rosetta low-quality and DMPFold models could be processed similarly.
```
find rosetta_high_quality_models -name "*.pdb.gz" | xargs -P 16 -n1 -I {} gunzip {}
```

## Running TM-align

TMAlign can be used to compare a single structure to a list of PDB formatted files. Running a single structure against all Rosetta or DMPFold models of MIP sequences should take under an hour.

The commands below assume that TMAlign has been installed and is in your `$PATH` and that `$(TARGET).pdb` is the PDB formatted model you want to search with. It first makes a list of models to operate on.

```
find rosetta_high_quality_models -name "*.pdb" | sed 's/rosetta_low_quality_models\///g' > pdb_list

TMalign ${TARGET}.pdb -dir2 rosetta_high_quality_models/ pdb_list -split 2 -ter 0 >& tmalign.log
```

## Processing TMAlign output

We have provided a script to parse the TMAlign output into a more convenient, tabular CSV format (`scripts/collect_tm_results.py`). The script takes the log file produced in the previous step and makes a file sorted based on the quality of the alignment. The most closely aligned structure will be at the end of the file.
```
python collect_tm_results.py -t tmalign.log -c tmalign.csv
```

# FUNCTION LOOKUPS OF DEEPFRI DATA

To search the DeepFRI function prediction data for MIP models we provide a script that allows for searching multiple functions at once and returning the results in a tabular CSV format given a DeepFRI score threshold (`scripts/function_to_MIP_id.py`). Score thresholds are between 0 and 1 with higher values indicating higher confidence. We typically use a score threshold of 0.1. The command line below will search for 10 functions in results for the rosetta high quality models with a threshold of 0.5 and use 16 processes at the same time. Each process requires about 1 GB of memory to store the precomputed DeepFRI results. As before the Rosetta low quality and DMPFold models could be processed similarly.

```
python function_to_MIP_id.py -t 0.50 -d rosetta_high_quality_function_predictions_test -n 16 -f GO:0043571 GO:0030683 GO:0004930 GO:0005125 GO:0009881 GO:0030246 EC:2.7.1.21 EC:2.7.10.- EC:4.2.1.1 EC:4.99.1.-
```

The command above produces the following files which contain MIP IDs, scores, GO/EC IDs, and GO/EC names

```
MIP_FUNCTIONS_BP_GO:0030683_0.5.csv
MIP_FUNCTIONS_BP_GO:0043571_0.5.csv
MIP_FUNCTIONS_EC_2.7.10.-_0.5.csv
MIP_FUNCTIONS_EC_2.7.1.21_0.5.csv
MIP_FUNCTIONS_EC_4.2.1.1_0.5.csv
MIP_FUNCTIONS_EC_4.99.1.-_0.5.csv
MIP_FUNCTIONS_MF_GO:0004930_0.5.csv
MIP_FUNCTIONS_MF_GO:0005125_0.5.csv
MIP_FUNCTIONS_MF_GO:0009881_0.5.csv
MIP_FUNCTIONS_MF_GO:0030246_0.5.csv
```
