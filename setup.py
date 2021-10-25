from setuptools import setup, find_packages
setup()

print("All set up! All scripts in the package can now be run from the command line with\noat <script_name>\nimporting functions in your own scripts can be done with\nfrom oxDNA_analysis_tools.<file name> import <function>")

print("")

print("If you would like to enable autocomplete in Bash, either add\n source oat-completion.sh\n to your .bashrc or run \nsudo cp oat-completion.sh /etc/bash_completion.d/\nto make autocompletes available to all users.")