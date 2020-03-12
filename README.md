# oxdna_analysis_tools

A suite of Python tools for performing generic structural analyses of oxDNA simulations.
Our goal in developing these tools is to provide a foundation of common analyses applicable to most simulations and to provide code examples for researchers looking to implement tools to meet their own research needs.
Running any script without arguments will print a brief description and list of required files.

An overarching description can be found in our preprint: https://doi.org/10.1101/2020.01.24.917419.


## Dependencies

This software package was written and tested using the following outside resources:<br/>
[Python](https://www.python.org/): 3.7 (minimum version 3.6),<br/>
[oxDNA](https://dna.physics.ox.ac.uk/index.php/Main_Page): 6985 (minimum version June 2019)<br/>
[NumPy](https://numpy.org/): 1.16,<br/>
[MatPlotLib](https://matplotlib.org/index.html): 3.0.3 (minimum version 3.0),<br/>
[BioPython](https://biopython.org/): 1.73,<br/>
[Scikit-Learn](https://scikit-learn.org/stable/): 0.21.2,<br/>
[Pathos](https://github.com/uqfoundation/pathos): 0.2.3</br>

To check if dependencies are correcly installed, run `config.py` with no arguments, this will check your environment and inform you which packages are incorrect.

## Brief script descriptions

Running instructions can be obtained for all scripts by running them with no arguments or the -h flag.

 * `align_trajectory.py` \<trajectory> \<topology> \<output> Aligns all configurations in the trajectory to the first configuration.<br/>
 * `backbone_flexibility.py` *WIP* (-p \<n_cpus>) \<trajectory> \<topology> \<output> Produces a color map of the deviation in backbone torsion angles over the course of the simulation.<br/>
 * `bond_analysis.py` (-p \<n_cpus>) \<input> \<trajectory> \<designed pairs file> \<output>  Produces and oxView overlay showing fractional occupancy of each designed pair.<br/>
 * `centroid.py` (-p \<n_cpus> -o \<centroid file> -i \<index file>) \<mean_structure> \<trajectory> \<topology> Takes a reference structure (usually a mean structure) and a trajectory and returns the structure in the trajectory with the lowest RMSF to the reference.<br/>
 * `clustering.py` \<serialized data input> Takes a set of configuration coordinates from other scripts and performs a DBSCAN clustering.  Produces trajectory files for each cluster and a visual representation of the clusters.  The -c option on pca.py and distance.py will call this script, which serializes its own data so you can re-launch the script to modify clustering parameters without re-running the analysis.<br/>
 * `compute_mean.py` (-p \<n_cpus> -f \<oxDNA/json/both> -o \<mean file> -d \<deviations file> -i \<index file> -a \<align conf id>) \<trajectory file> \<topology> Produces the mean structure of a trajectory via single-value decomposition superposition.  If the -i flag is added with a file containing a list of particle IDs, the mean structure will be calculated based only on the subset of particles included in the list.  These lists can be downloaded from oxView using the "Download selected base list" button.  By default, this script aligns to a random configuration in the trajectory.  However, if you would like to align to a specific configuration, you can specify its position in the trajectory with the -a flag.  The -d flag will automatically run compute_deviations.py from the mean structure.<br\>
 * `compute_deviations.py`(-p \<n_cpus> -o \<deviations file>) \<mean structure> \<trajectory file> \<topology> Computes the per-nucleotide RMSF from the mean structure.  Can be called automatically by compute_mean with the -d option.
 * `config.py` Contains system specific information.  Update before you start running these scripts.  Running with no arguments will run a dependency check on your environment.
 * `contact_map.py` (-v) \<input> \<trajectory> produces a visual contact map of internucleotide distances if the -v option is given.  Otherwise lists all distances.<br/>
 * `distance.py` (-o \<output> -f \<histogram/trajectory/both> -d \<data file output>) -i \<\<input> \<trajectory> \<particleID 1> \<particleID 2>> The -i option can be called as many times as needed.  Additional calls will be overlaid on the same graph. Produces the user's choice of histograms, timeseries or text outputs (set by -f and -d options).<br/>
 * `duplex_angle_finder.py` (-p \<n_cpus> -o \<output>) \<input> \<trajectory> Produces a text file containing identification information for all duplexes at each configuration in the trajectory.<br/>
 * `duplex_angle_plotter.py` (-o \<output> -f \<histogram/trajectory/both>) -i \<\<angle file> \<particleID 1> \<particleID 2>> Produces either histograms or timeseries of the angle between specified duplexes.<br/>
 *  `eRMSD.py` (-p \<n_cpus>) \<input> \<trajectory> Computes the [eRMSD](https://academic.oup.com/nar/article/42/21/13306/2903225) between each pair of configurations.  This is a very computationally intense operation and not recommended for large structures or large trajectories. <br/>
 * `forces2pairs.py` \<force file> (\<output>) Takes an oxDNA external forces file and produces a list of space-separated pairs.  This output is used as an input for `bond_analysis.py`. <br/>
 * `generate_force.py` (-o \<output> -f \<pairs file>) \<input> \<configuration> Produces an external force file enforcing the current base pair arrangement.  Useful for generating pairs files using `forces2pairs.py` or for enforcing relaxation to a particular configuration during an oxDNA simulation.<br/>
 * `mean2dat.py` \<mean.json> \<output>  Converts the .json configuration file produced by `compute_mean.py` to a .dat file readable by oxDNA or oxView.<br/>
 * `multidimensional_scaling_mean.py` (-p \<n_cpus>) \<input> \<trajectory> \<mean prefix> \<dev prefix> Computes the mean structure based on local pairwise distances An alternative to `compute_mean.py` that works better for highly flexible structures.<br/>
 * `output_bonds.py` (-v \<output>) \<input> \<trajectory> Lists all the interactions between nucleotides.  The output is the same as the `output_bonds.py` found in the oxDNA UTILS, this one is written in Python3, and has the option of an oxView overlay showing average total energy per nucleotide. <br/>
 * `pca.py` (-p \<n_cpus> -c) \<input> \<trajectory> \<mean file> \<output> Computes the principal components of deviations away from the mean.  The principal components are written as an oxView overlay file.  If the -c flag is used it will also run the clustering script on each configuration's position in principal component space. <br/>
 * `superimpose.py` \<topology> \<configuration> \<configuration> (<configuration> <configuration> ...) Superimposes all further configurations onto the first file provided.
 
### UTILS
The UTILS directory contains utility modules used by other scripts in this package.

* `base.py` A python3 update of the `base.py` script found in the oxDNA distribution.  This contains class definitions for nucleotide/strand/system.  These are used to create, modify and write oxDNA systems in a Python environment. <br/>
* `geom.py` A set of algorithms to find various geometric parameters of DNA/RNA helices.  Currently only the axis fitting function is used. <br/>
* `model.h` The model parameters of the oxDNA model.  Used by base.py.
* `parallelize.py` The parallelization module used by the analysis scripts.  Splits the trajectory file into temporary files and attaches a reader managed by a different CPU to each chunk. <br/>
* `parallelize_old.py` An older implementation of the parallelizer that does not split into temporary files.  Certain multiprocessing architectures can have issues with reader congestion when using this scheme.  Left in place in case of use cases where memory usage is an issue.
* `readers.py` Contains utility functions for working with oxDNA files, including extracting input file parameters, calculating the number of configurations in a trajectory and creating a system as defined in `base.py` from a configuration/topology pair.

## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

## A note on defaults
The config.py script needs to be modified to contain your own path to DNAnalysis.  If you are pushing modifications back to GitHub, you need to set `git update-index --skip-worktree config.py` so that your changes don't overwrite other user's path when they do a pull. See:  
https://stackoverflow.com/questions/4348590/how-can-i-make-git-ignore-future-revisions-to-a-file

The function used to bring structures back in the simulation box, `inbox_system()` in UTILS/base.py requires a reference nucleotide to overcome the combination of periodic boundary conditions and fix_diffusion artifacts.  This nucleotide should be something in the center of the structure as a particle on the edge can cause improper inboxing.  The ID of the reference is set in config.py, so if you encounter issues with your structure and periodic boundaries, change this number.

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
Erik Poppleton, Joakim Bohlin, Michael Matthies, Shuchi Sharma, Fei Zhang, Petr Sulc: Design, optimization, and analysis of large DNA and RNA nanostructures through interactive visualization, editing, and molecular simulation (bioarxiv preprint: https://doi.org/10.1101/2020.01.24.917419)
