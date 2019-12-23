# oxdna_analysis_tools

A suite of Python tools for performing generic structural analyses of oxDNA simulations.
Our goal in developing these tools is to provide a foundation of common analyses applicable to most simulations and to provide code examples for researchers looking to implement tools to meet their own research needs.
Running any script without arguments will print a brief description and list of required files.

More detailed descriptions can be found at INSERT PAPER LINK HERE.


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

 * `contact_map.py` Takes an input and configuration file and produces a visual contact map of internucleotide distances<br/>
 * `align_trajectory.py` Takes a trajectory and topology file pair and aligns all configurations in the trajectory to the first configuration.<br/>
 * `all_vectors.py` Takes an input and trajectory file pair and produces a matrix of vectors between every nucleotide at each step.<br/>
 * `backbone_flexibility.py` *WIP* Takes an input and trajectory file pair and produces a color map of the deviation in backbone torsion angles over the course of the simulation.<br/>
 * `bond_analysis.py` Takes and input and trajectory file pair along with a text file specifying the designed bonded pairs.  Produces and oxView overlay showing fractional occupancy of each designed pair.<br/>
 * `centroid.py` *WIP* Takes a set of configuration coordinates from other scripts and produces the configuration file for the centroid along that coordinate.<br/>
 * `clustering.py` Takes a set of configuration coordinates from other scripts and performs a DBSCAN clustering.  Produces trajectory files for each cluster and a visual representation of the clusters.<br/>
 * `distance.py` Takes any number of input and trajectory file pairs along with specific nucleotide pairs to find distances between and produces the user's choice of histograms, timeseries or text outputs.<br/>
 * `duplex_angle_finder.py` Takes an input and trajectory file pair and produces a text file containing identification information for all duplexes at each configuration in the trajectory.<br/>
 * `duplex_angle_plotter.py` Takes the duplex text file from `duplex_angle_finder.py` and produces either histograms or timeseries of the angle between specified duplexes.<br/>
 *  

## Output files and visualization

Many scripts in this package produce data overlay json files that can be used with [oxView](https://github.com/sulcgroup/oxdna-viewer).
To load an overlay, drag and drop the json file along with the configuration and topology files, or drag and drop the json file once the load is completed.

## Citation

If you use these scripts or oxView in your published work, please cite:<br/>
PAPER CITATION HERE
