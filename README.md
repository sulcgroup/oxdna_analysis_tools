# oxdna_analysis_tools

A suite of Python tools for performing generic structural analyses of oxDNA simulations.
Our goal in developing these tools is to provide a foundation of common analyses applicable to most simulations and to provide code examples for researchers looking to implement tools to meet their own research needs.
Running any script without arguments will print a brief description and list of required files.

An overarching descript can be found at INSERT PAPER LINK HERE.


## Dependencies

This software package was written and tested using the following outside resources:<br/>
[Python](https://www.python.org/): 3.7 (minimum version 3.6),<br/>
[oxDNA](https://dna.physics.ox.ac.uk/index.php/Main_Page): 6985 (minimum version June 2019)<br/>
[NumPy](https://numpy.org/): 1.16,<br/>
[MatPlotLib](https://matplotlib.org/index.html): 3.0.3 (minimum version 3.0),<br/>
[BioPython](https://biopython.org/): 1.73,<br/>
[Scikit-Learn](https://scikit-learn.org/stable/): 0.21.2,<br/>
[Pathos](https://github.com/uqfoundation/pathos): 0.2.3

## Brief script descriptions

Running instructions can be obtained for all scripts by running them with no arguments or the -h flag.

 * `contact_map.py` Takes an input and configuration file and produces a visual contact map of internucleotide distances<br/>
 * `align_trajectory.py` Takes a trajectory and topology file pair and aligns all configurations in the trajectory to the first configuration.<br/>
 * `all_vectors.py` Takes an input and trajectory file pair and produces a matrix of vectors between every nucleotide at each step.<br/>
 * `backbone_flexibility.py` *WIP* Takes an input and trajectory file pair and produces a color map of the deviation in backbone torsion angles over the course of the simulation.<br/>
 * `bond_analysis.py` Takes an input and trajectory file pair along with a text file specifying the designed bonded pairs.  Produces and oxView overlay showing fractional occupancy of each designed pair.<br/>
 * `centroid.py` *WIP* Takes a set of configuration coordinates from other scripts and produces the configuration file for the centroid along that coordinate.<br/>
 * `clustering.py` Takes a set of configuration coordinates from other scripts and performs a DBSCAN clustering.  Produces trajectory files for each cluster and a visual representation of the clusters.<br/>
 * `distance.py` Takes any number of input and trajectory file pairs along with specific nucleotide pairs to find distances between and produces the user's choice of histograms, timeseries or text outputs.<br/>
 * `duplex_angle_finder.py` Takes an input and trajectory file pair and produces a text file containing identification information for all duplexes at each configuration in the trajectory.<br/>
 * `duplex_angle_plotter.py` Takes the duplex text file from `duplex_angle_finder.py` and produces either histograms or timeseries of the angle between specified duplexes.<br/>
 *  `eRMSD.py` Takes an input and trajectory file and computes the [eRMSD](https://academic.oup.com/nar/article/42/21/13306/2903225) between each pair of configurations.  This is a very computationally intense operation and not recommended for large structures or large trajectories. <br/>
 * `forces2pairs.py` Takes an oxDNA external forces file and produces a list of space-separated pairs.  This output is used as an input for `bond_analysis.py`. <br/>
 * `generate_force.py` Takes an input and configuration file and produces an external force file enforcing the current base pair arrangement.  Useful for generating pairs files using `forces2pairs.py` or for enforcing relaxation to a particular configuration during an oxDNA simulation.<br/>
 * `mean2dat.py` Converts the .json configuration file produced by `compute_mean.py` to a .dat file readable by oxDNA or oxView.<br/>
 * `multidimensional_scaling_mean.py` Takes an input and configuration file and computes the mean structure based on local pairwise distances An alternative to `compute_mean.py` that works better for highly flexible structures.<br.>
 * `output_bonds.py` Takes an input and configuration file and lists all the interactions between nucleotides.  The output is the same as the `output_bonds.py` found in the oxDNA UTILS, this one is written in Python3, has a parallelization option and has the option of an oxView overlay showing average total energy per nucleotide. <br/>
 * `pca.py` Takes an input, a mean structure, and a configuration file and computes the principal components of deviations away from the mean.  The principal components are written as an oxView overlay file.  Also has the option to run the clustering script on each configuration's position in principal component space. <br/>
 * `superimpose.py` Takes a topology file and any number of configuration files that correspond to that topology and superimposes all further configurations onto the first file provided.
 
### UTILS
The UTILS directory contains utility modules used by other scripts in this package.

* `base.py` A python3 update of the `base.py` script found in the oxDNA distribution.  This contains class definitions for nucleotide/strand/system.  These are used to create, modify and write oxDNA systems in a Python environment. <br/>
* `geom.py` A set of algorithms to find various geometric parameters of DNA/RNA helices.  Currently only the axis fitting function is used. <br/>
* `model.h` The model parameters of the oxDNA model.  Used by base.py.
* `parallelize.py` The parallelization module used by the analysis scripts.  Splits the trajectory file between a given number of workers and runs the provided function on each chunk of the trajectory file. <br/>
* `readers.py` Contains utility functions for working with oxDNA files, including extracting input file parameters, calculating the number of configurations in a trajectory and creating a system as defined in `base.py` from a configuration/topology pair.

## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
PAPER CITATION HERE
