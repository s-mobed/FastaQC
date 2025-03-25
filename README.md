This app creates two versions of FASTA file based of the output FASTA file from Piranha. 
The first is the QC_PASS version which contains only consensus FASTAs from samples that passed QC and will later be reported to the Polio Consortium.
The second is Phylogenetic Database version, where the header is modified to give clearer information and contains all samples from run excluding controls. This will be used to populate the phylo trees in the Piranha Reports, when using the phylo module in Piranha. Allowing the user to indentify possible contamination with prior samples and associate Vaccine Derived Polio viruses (VDPVs) to certain known lineages.
