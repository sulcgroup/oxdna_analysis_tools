# oxDNA Analysis Tools

A suite of Python tools for performing generic structural analyses of oxDNA simulations.
Our goal in developing these tools is to provide a foundation of common analyses applicable to most simulations and to provide code examples for researchers looking to implement tools to meet their own research needs.
Running any script without arguments will print a brief description and list of required files.

An overarching description can be found in this paper: https://academic.oup.com/nar/article/48/12/e72/5843822.


## Dependencies and installation

### Pip installation
oxDNA analysis tools can be installed from PyPi via pip:

`pip install oxDNA-analysis-tools`

This will also install all dependencies.  Bash autocompletes will not be set up, see below for setting up autocompletion.

### Installation from source
It can also be installed from the [GitHub repository](https://github.com/sulcgroup/oxdna_analysis_tools) or the zip file of the source code available on PyPi via the following method:  

1. Clone the repository or download and inflate the zip file.  
2. Run one of the following commands (pip to automatically install dependencies or setup.py if you would like to manage them yourself):  
   `pip install .`  
   `python setup.py install`  

If you are not installing via pip, the following dependencies are required and can all be obtained from either pip or conda:  
[Python](https://www.python.org/): 3.7 (minimum version 3.6),<br/>
[oxDNA](https://dna.physics.ox.ac.uk/index.php/Main_Page): 6985 (minimum version June 2019)<br/>
[NumPy](https://numpy.org/): 1.16,<br/>
[MatPlotLib](https://matplotlib.org/index.html): 3.0.3 (minimum version 3.0),<br/>
[BioPython](https://biopython.org/): 1.73,<br/>
[Scikit-Learn](https://scikit-learn.org/stable/): 0.21.2,<br/>
[Pathos](https://github.com/uqfoundation/pathos): 0.2.3</br>

### Setting up Bash autocompletes
The invocation `oat` is calling a Python script which then handles calling the other available scripts.  If you would like autocompletes for the specific script names (and are using a Unix command line), these are provided by `oat-completion.sh` which can also be found in the repository.  To add autocompletes to your system, either append it to your local `.bash_completion` file with:  

`cat oat-completion.sh >> ~/.bash_completion`

Or add it to your global completions with:  

`sudo cp oat-completion.sh /etc/bash_completion.d/`

### DNAnalysis and verifying installation
Some scripts require DNAnalysis, a program in the oxDNA distribution to calcuate the oxDNA energy function for identifying bonded nucleotides.  Unfortunatley, there is not a good way to detect the location of this program on each user's system when installing from pip, so the path to it must be hardcoded.  In order to use scripts which check base pairing (particularly `bond_analysis`, `multidimensional_scaling_mean`, and `duplex_angle_finder`), you need to locate `config.py` in the oxDNA_analysis_tools source directory. There, edit the `PROCESSPROGRAM` variable in to point to the compiled DNAnalysis binary.  

To verify your environment is set up correctly, run 

`oat config`

The scripts which use DNAnalysis require at least release 6985 (June 2019) of oxDNA. Unfortunately, there is not currently a good way to check oxDNA version (this has been fixed in the [bleeding edge repository](https://github.com/lorenzo-rovigatti/oxDNA)).  If `multidimensional_scaling_mean.py` or `duplex_angle_finder.py` throw an error from oxDNA, then your version is probably out of date.

## Using oxDNA analysis tools
Once installed, all standalone scripts can be called from the command line via the following invocation:  
`oat <script name> <script arguments>`  

For example, to compute the mean structure and deviations of a file called `trajectory.dat` using 4 CPUs and outputting to files called `mean.dat` and `devs.json`, you would run:  
`oat compute_mean -p 4 -o mean.dat -d devs.json trajectory.dat`

To see a detailed description of the script command line arguments, run the script with the `-h` flag.

These scripts are intended to be extensible and re-used for custom analysis by users.  The functions in this library can be imported into your Python scripts via:  
`from oxDNA_analysis_tools.<script name> import <object name>`

So for example, if you would like to use the ErikReader trajectory reader, you would include this line at the start of your Python script:  
`from oxDNA_analysis_tools.UTILS.readers import ErikReader`

## File formats

This package mostly uses the oxDNA files as described in [the oxDNA documentation](https://dna.physics.ox.ac.uk/index.php/Documentation).  A brief description of each file is provided here for easy reference:  
**trajectory** - A file containing a sequence of oxDNA configurations.  Each configuration starts with a three line header containing the timestep, box size and energy information.  There is then one line per particle with 15 values corresponding to the position, orientation and velocity of each particle.  
**topology** - A file containing sequence and connectivity information of the simulated structure.  The first line defines the number of particles and the number of strands.  There is then one line per particle with 4 values corresponding to the strand ID, the base type, the 3' neighbor and the 5' neighbor of each particle.  Note that oxDNA files are written 3'-5' rather than the traditional 5'-3'.  
**input** - The input file used to run oxDNA.  This contains simulation information such as number of steps, simulation method and temperature as well as I/O information.  Example files can be found in the "example_input_files" and "paper_examples" directories.  
**force file**: An oxDNA mutual trap file that defines an external force between two particles in a simulation.  This is also defined in [the oxDNA documentation](https://dna.physics.ox.ac.uk/index.php/Documentation).  

The Following files are unique to this package:  
**oxView json file**: This file contains overlay information that can be read by [oxView](https://github.com/sulcgroup/oxdna-viewer).  There are two different formats that are produced by these scripts.  The first is a one-value-per-particle file that creates a colormap overlay with extreme colors corresponding to the minimum and maximum values in the file.  The second is a three-values-per-particle file that oxView uses to draw arrows on the scene.  OxView automatically differentiates files based on the number of values corresponding to each particle.  
**designed pairs file**: This file contains a list of particle pairs in the intended design.  Each line corresponds to a pair and each pair is a space-separated pair of particle IDs.  Designed pairs files can be generated by forces2pairs.py and generate_force.py.  
**angle file**: The output file generated by duplex_angle_finder.py.  Details on the format can be found in a comment in the duplex_angle_plotter.py script, but briefly each line contains starting and ending nucleotides and orientation data for a duplex in the structure.  Like trajectories, this contains information for every configuration in a trajectory.  
**index file**: A space-seperated list of particle IDs that can be used by compute_mean.py for subset alignment.  It can be generated by the "Download Selected Base List" button in oxView.  
**serialized data input**: To make it easy to adjust clustering parameters, the clustering script serializes its input in json format so the script can be re-launched quickly with this file as the only input.  


## Brief script descriptions

Running instructions can be obtained for all scripts by running them with no arguments or the -h flag.

 * `align_trajectory.py` (-i \<index_fine> -r\<configuration_file>) \<trajectory> \<output> Aligns all configurations in the trajectory to the first configuration.  Produces an oxDNA trajectory that is a copy of the provided trajectory with all translations and rotations removed.  If the -i flag is added with an index file containing a list of particle IDs, the alignment will be calculated based only on the subset of particles included in the list.  These lists can be created from selected bases in oxView using the "Download selected base list" button.  If the -r flag is added with a configuration file, the trajectory will be aligned to the provided configuration rather than the first configuration in the trajectory.<br/>
 * `backbone_flexibility.py` *WIP* (-p \<n_cpus>) \<trajectory> \<topology> \<output> Produces an oxView json color map of the deviation in backbone torsion angles over the course of the simulation.<br/>
 * `bond_analysis.py` (-p \<n_cpus>) \<input> \<trajectory> \<designed pairs file> \<output>  Calculates the hydrogen bond occupancy compared with the intended design.  Produces an oxView json file to \<output> that creates a visual colormap corresponding to the occupancy of each pair observed in the simulation.<br/>
 * `centroid.py` (-p \<n_cpus> -o \<centroid file> -i \<index file>) \<mean_structure> \<trajectory> Takes a reference structure (usually a mean structure) and a trajectory and returns the structure in the trajectory with the lowest RMSF to the reference as an oxDNA configuration file.<br/>
 * `clustering.py` \<serialized data input> Takes a set of configuration coordinates from other scripts and performs a DBSCAN clustering.  Produces trajectory files for each cluster (note that these trajectories do not necessarily contain contiguous timesteps) and a visual representation of the clusters in either a 2D or 3D plot (2D if only 2 dimensions of input data are given, 3D if 3 or more are given with the first three displayed).  The -c option on pca.py and distance.py will call this script. Clustering.py serializes its own data to a file called `cluster_data.json` so you can re-launch the script to modify clustering parameters without re-running the analysis.  Clustering parameters (`EPS` and `MIN_SAMPLES` can be modified in the script at the start of the `perform_DBSCAN()` function.<br/>
 * `compute_mean.py` (-p \<n_cpus> -f \<oxDNA/json/both> -o \<mean structure> -d \<deviations file> -i \<index file> -a \<align conf id>) \<trajectory file> Produces the mean structure as either an oxDNA configuration file or as a json file that contains the structure broken down into positions and rotations of a trajectory via single-value decomposition superposition.  If the -i flag is added with an index file containing a list of particle IDs, the mean structure will be calculated based only on the subset of particles included in the list.  These lists can be created from selected bases in oxView using the "Download selected base list" button.  By default, this script aligns to a random configuration in the trajectory.  However, if you would like to align to a specific configuration, you can specify its position in the trajectory with the -a flag.  The -d flag will automatically run compute_deviations.py from the mean structure.<br/>
 * `compute_deviations.py`(-i\<index_file> -p \<n_cpus> -o \<deviations file> -r \<rmds plot> -d \<order parameter file>) \<mean structure> \<trajectory file> Computes the per-nucleotide RMSF from the mean structure.  Can be called automatically by compute_mean with the -d option. Produces an oxView json flie that colors each particle based on its RMSF. If the -i flag is added with an index file containing a list of particle IDs, the deviations will be calculated based only on aligning the subset of particles included in the list.  These lists can be created from selected bases in oxView using the "Download selected base list" button. If the -r option is called with a plot file name, the deviations over time will be displayed as a line plot.  If the -d option is called with an output file name, the script will produce an oxView order parameter file which can be loaded in the "Trajectory" tab.<br/>
 * `config.py` Contains system specific information and performs dependency checks.  Update the path to your compiled DNAnalysis in this script before you start running these scripts.  Running with no arguments will run a dependency check on your environment, which is recommended after downloading.
 * `contact_map.py` (-v) \<input> \<trajectory> produces a contact map of internucleotide distances if the -v option is given.  Otherwise lists all distances.<br/>
 * `distance.py` (-c -o \<output> -f \<histogram/trajectory/both> -d \<data file output>) -i \<\<input> \<trajectory> \<particleID 1> \<particleID 2> (\<particleID 1> \<particleID 2> ...)> Computes the distance between provided particle pairs. The -i option can be called multiple times to overlay data from multiple trajectories.  Additional calls will be overlaid on the same graph. Produces the user's choice of histograms, timeseries or text outputs (set by -f and -d options).  Data series names can be set by modifying the `names` variable in the script.  The -c option will run the output of the distances through the clusterin script<br/>
 * `duplex_angle_finder.py` (-p \<n_cpus> -o \<output>) \<input> \<trajectory> Produces an angle file containing identification information for all duplexes at each configuration in the trajectory.  This file is visualized by duplex_angle_plotter.py<br/>
 * `duplex_angle_plotter.py` (-o \<output> -f \<histogram/trajectory/both>) -i \<\<angle file> \<particleID 1> \<particleID 2>> Reads the angle file produced by duplex_angle_finder and produces either histograms or timeseries of the angle between specified duplexes.  Can provide the same angle file multiple times with different particle IDs or different angle files by calling -i multiple times.  The data from all calls will be drawn on the same graph. Data series names can be set by modifying the `names` variable in the script<br/>
 *  `eRMSD.py` (-p \<n_cpus>) \<input> \<trajectory> Computes the [eRMSD](https://academic.oup.com/nar/article/42/21/13306/2903225) between each pair of configurations.  This is a very computationally intense operation and not recommended for large structures or large trajectories. <br/>
 * `forces2pairs.py` \<force file> (\<output>) Takes an oxDNA external forces file and produces a pairs file.  This output is used as an input for `bond_analysis.py`. <br/>
 * `generate_force.py` (-o \<output> -f \<pairs file>) \<input> \<configuration> Produces an external force file enforcing the current base pair arrangement. Useful for generating pairs files using `forces2pairs.py` or for enforcing relaxation to a particular configuration during an oxDNA simulation.  The -f option will automatically create the pairs file.<br/>
 * `mean2dat.py` \<mean.json> \<output>  Converts the .json configuration file produced by `compute_mean.py` to a .dat file readable by oxDNA or oxView.<br/>
 * `multidimensional_scaling_mean.py` (-p \<n_cpus>) \<input> \<trajectory> \<mean prefix> \<dev prefix> Computes the mean structure based on local pairwise distances between particles.  An alternative to `compute_mean.py` that works better for highly flexible structures.  Produces an oxDNA configuration/topology file pair and an oxView json file showing the per-particle deviation in distance to neighboring particles.<br/>
 * `output_bonds.py` (-v \<output>) \<input> \<trajectory> Lists all the interactions between nucleotides.  The output is the same as the `output_bonds.py` found in the oxDNA UTILS, this one is written in Python3, and has the option of an oxView json overlay showing average total energy per nucleotide. <br/>
 * `oxDNA_PDB.py` (topology configuration_file direction -p <list of protein pdb files in system in order of occurence>)  Converts either oxDNA DNA files to pdb, or oxDNA DNA/Protein hybrids to pdb format.br/>
 * `pca.py` (-p \<n_cpus> -c) \<input> \<trajectory> \<mean file> \<output> Computes the principal components of deviations away from the mean.  The principal components are written as an oxView json overlay file that will show the direction of the top mode.  More components can be summed and added to the overlay by modifying the "SUM" variable in the script.  If the -c flag is used it will also run the clustering script on each configuration's position in principal component space. <br/>
 * `superimpose.py` (-i \<index_file>) \<configuration> \<configuration> (<configuration> <configuration> ...) Superimposes all further configurations onto the first file provided.  Produces an oxDNA trajectory file that is a copy of the input file with all translations and rotations removed.  It is expected that all referenced nucleotides are the same, so the configurations either must share the same topology, or you must only align to the shared particles using an index file. An index file can be downloaded from oxView using the "Download selected base list" button.
 
### UTILS
The UTILS directory contains utility modules used by other scripts in this package.

* `all_vectors.py` A wrapper for the all_vectors observable in oxDNA that computes the relative position of all particles in a configuration. <br/>
* `base.py` A python3 update of the `base.py` script found in the oxDNA distribution.  This contains class definitions for nucleotide/strand/system.  These are used to create, modify and write oxDNA systems in a Python environment. <br/>
* `geom.py` A set of algorithms to find various geometric parameters of DNA/RNA helices.  Currently only the axis fitting function is used. <br/>
* `model.h` The model parameters of the oxDNA model.  Used by base.py. <br/>
* `parallelize.py` The parallelization module used by the analysis scripts.  Splits the trajectory file into temporary files and attaches a reader managed by a different CPU to each chunk.  Each reader then feeds successive configurations into a specified function.<br/>
* `parallelize_old.py` An older implementation of the parallelizer that does not split into temporary files.  Certain multiprocessing architectures can have issues with reader congestion when using this scheme.  Left in place in case of use cases where memory usage is an issue.
* `pdb.py` Helper Functions/Classes for pdb conversion <br/>
* `protein_to_pdb` Contains protein specific functions for protein to pdb conversion<br/>
* `readers.py` Contains utility functions for working with oxDNA files, including extracting input file parameters, calculating the number of configurations in a trajectory and creating a system as defined in `base.py` from a configuration/topology pair.<br/>
* `utils.py` Contains utility functions for pdb conversion<br/>
* `dd12_na.pdb` Used during pdb conversion script
## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

By default scripts in this package that produce graphs save them as .png files.  All graphing is done using the Matplotlib interface and users are encouraged to make modifications to the graph styles to fit their unique needs.

## A note on defaults
The config.py script needs to be modified to contain your own path to DNAnalysis.  If you are pushing modifications back to GitHub, you need to set `git update-index --skip-worktree config.py` so that your changes don't overwrite other user's path when they do a pull. See:  
https://stackoverflow.com/questions/4348590/how-can-i-make-git-ignore-future-revisions-to-a-file

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
Erik Poppleton, Joakim Bohlin, Michael Matthies, Shuchi Sharma, Fei Zhang, Petr Sulc: Design, optimization, and analysis of large DNA and RNA nanostructures through interactive visualization, editing, and molecular simulation. (2020) Nucleic Acids Research e72. https://doi.org/10.1093/nar/gkaa417
