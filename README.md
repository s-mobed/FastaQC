This app creates two versions of FASTA file based of the output FASTA file from Piranha. 
The first is the QC_PASS version which contains only consensus FASTAs from samples that passed QC and will later be reported to the Polio Consortium.
The second is Phylogenetic Database version, where the header is modified to give clearer information and contains all samples from run excluding controls. This will be used to populate the phylo trees in the Piranha Reports, when using the phylo module in Piranha. Allowing the user to indentify possible contamination with prior samples and associate Vaccine Derived Polio viruses (VDPVs) to certain known lineages.

# Installation

Download the compiled version found in the releases section (https://github.com/s-mobed/FastaQC/releases).

# Usage

<p align="center" width="100%">
    <img width="75%" src="https://github.com/user-attachments/assets/a4f6d08e-0c4d-4333-b138-70d6684f5741">
</p>

1. The Dropbox is used as the input for one detailed run report CSV file and its associated FASTA file. It will not accept more than two files and file extensions
   other than .csv and .fasta
2. The Destination text box will show the chosen destination of the output FASTA files.
3. These four buttons are:
   a. Destination button will bring a folder selection menu that allows the user to choose the destination of the files. This will be displayed above.
   b. Parse files button starts the underlying code of the app to check the detailed run report and FASTA file to see which samples have results and which have
      passed or failed. It will also highlight any inconsistencies between them.
   c. Generate QC Fasta button will be grayed out until the files have been parsed. Once parsed and the user has checked the displayed results, the user can press
      to output the files.
   d. Clear button will clear the selected files to allow the user to process further pairs of samples.

4. The Method used for the run will be written in the header of the DB FASTA versions. By default it is DDNS, but the user can choose between both methods. 
