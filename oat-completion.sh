#/usr/bin/env bash

_oat_completions(){
	if [ "${#COMP_WORDS[@]}" != "2" ]; then
		_longopt
		return
	fi


	COMPREPLY=($(compgen -W "align_trajectory anm_parametrize backbone_flexibility bond_analysis centroid clustering compute_deviations compute_mean config contact_map distance duplex_angle_finder duplex_angle_plotter eRMSD forces2pairs generate_force mean2dat minify multidimensional_scaling_mean output_bonds p_align pca pca_experimental plot_energy subset_trajectory superimpose" "${COMP_WORDS[COMP_CWORD]}"))
}

complete -F _oat_completions oat
