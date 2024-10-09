# Project Pair Master analysis

Much of this analysis is available in VDJServer's RepCalc workflow, but these
scripts were written before that was available. At some point, we should replace
these scripts with output from RepCalc on VDJServer. Right now this workflow is
manually running RepCalc.

## Analysis Workflow

When the OrderNumber_Original folder with OrderNumber_TSVready.xls in OneDrive is ready for analysis, the following workflow describes the steps involved to generate annotated sequences. This is primarily technical information to remind Scott what programs and tools to use so the same reproducible analysis is performed.

1. Copy folder with OrderNumber_TSVready.xls workbook to local drive (data/immune/MonsonLab/Sanger). Avoid complications with OneDrive and prevent polluting with intermediate and/or temporary analysis files.
2. Export OrderNumber_TSVready.xls to CSV file (not UTF-8 CSV) using Excel. Might need to change the names of files and/or folders in case they contain spaces, just to make them easier to process at command line. Columns are “PCR product name” and “Sequence”.
3. Use the docker-repcalc alias to start docker container for running local programs.
   1. cd ~/Projects
   2. docker-repcalc
   3. cd /work/data/immune/MonsonLab/Sanger/FOLDER
4. Generate FASTA file from the exported CSV. Provide the script with the CSV filename and the output filename. Produce NAME.fasta where NAME is the name of the folder for the data.
  a. python3 /work/research/immune/code/MonsonLab/RHAB/make_fasta.py CSV NAME
  b. Verify the number of sequences in the FASTA file and that they all have unique sequence IDs
5. Upload FASTA file to VDJServer and run IgBlast (Human/IG). I've been putting the files in a single project named: Monson Lab B cells (Sanger SC). It is useful to name the job with NAME to easily match the results later.
6. Download the Archive Output from IgBlast and unzip files in the vdjserver folder on local drive. The name of the folder with the unzipped output is the job ID. The AIRR TSV files within the archive are gzipped, so ungzip them. We are mainly interested in *.makedb.airr.tsv and *.igblast.airr.tsv
7. Generate final analysis. The final analysis steps have been combined into a single script.
  a. bash /work/research/immune/code/MonsonLab/RHAB/rhab_summary_v4.sh NAME
  b. This script performs the operations:
    i. Calculate CDR3 properties.
    ii. Create germline alignments with IMGT gaps and D mask.
    iii. Calculate mutations.
    iv. Generate summary TSV.
8. Copy files back to OneDrive.
  a. FASTA: seqs-FOLDER.fasta
  b. SUMMARY: seqs-FOLDER.summary.tsv
  c. IgBlast AIRR TSV: seqs-FOLDER.igblast.airr.tsv
  d. Mutations AIRR TSV: seqs-FOLDER.mutations.airr.tsv
  e. MonsonLab_$NAME.mutations.airr.tsv
