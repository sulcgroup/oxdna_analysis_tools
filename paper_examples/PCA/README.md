This directory contains an example to compute the principal components of a holliday junction as shown in figure 10 of the paper.

1. Run an oxDNA simulation using the provided input file.

   /path/to/oxDNA input_dna

   A couple of notes on this simulation:
     This simulation is small and can be run on a laptop or non-high end computer.
     The input file here will run for 1e8 steps.  For production runs we usually recomend 1e9 steps, however 1e8 will be fine for an example.
     This simulation uses the sequence-dependent model of DNA.  By default, oxDNA uses an average sequence model, but in cases where the different strengths of interactions between the four nucleobases may be important, it is more appropriate to use sequence dependence (see Snodin et. al. (2015)).

2. Compute the mean structure via SVD

   python ../../compute_mean.py -f oxDNA -o mean.dat trajectory_trap.dat holliday.top

   PCA is defined by deviations from some reference.  In this case, we will use the mean structure as a reference.

3. Compute the principal components

   python ../../pca.py input_dna trajectory_trap.dat mean.dat pca.json

   A few notes on the PCA script:
   Because of the fix diffusion process that oxDNA uses to maintain correct physics, ensuring that a structure is not cut by the periodic boundary is a non-trivial task.  To solve this problem, some analysis scripts require a reference nucleotide which everything else is centered around. By default, this is nucleotide 0.  In this case, nucleotide 0 is fine, however in other structures, if you observe an area with very high RMSF in the mean structure or unreasonably high displacement in the PCA, change the INBOXING_REFERENCE_PARTICLE in config.py.
   If you are running on a computer with multiple CPUs, you can additionally specify the -p <number> option to compute principal components in parallel using tht many CPUs.  This results in significant performance gains up to around 30 cpus
   The -c option will run the clustering algorithm on reconstructions of each configuration in principal component space.  By default this uses all components, but if you want to use only the top few components, uncomment the line that truncates the linear terms in the last block of code in the script.

This will produce an oxView json file that plots arrows on the structure which corresond to the weighted sum of the first n components.  n can be set by modifying the SUM variable in the pca.py script.  To view in the viewer, drag and drop the mean structure, topology and PCA json files onto the viewer window. 

This script also creates a scree plot showing the weights of components.