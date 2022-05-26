#/usr/bin/env bash

_oat_completions(){
	if [ "${#COMP_WORDS[@]}" != "2" ]; then
		_longopt
		return
	fi


	COMPREPLY=($(compgen -W "align anm_parameterize backbone_flexibility bond_analysis centroid clustering config contact_map db_to_force deviations distance duplex_angle_plotter duplex_finder forces2pairs generate_force mean minify multidimensional_scaling_mean output_bonds oxDNA_PDB pca plot_energy subset_trajectory superimpose" "${COMP_WORDS[COMP_CWORD]}"))
}

complete -F _oat_completions oat
