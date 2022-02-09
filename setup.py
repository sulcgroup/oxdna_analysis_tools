from setuptools import setup, find_packages
import re
import os

def get_property(prop):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop), open('src/oxDNA_analysis_tools/__init__.py').read())
    return result.group(1)

setup(version = get_property('__version__'))

print("\n\n################################################################################")
print("\nAll set up! All scripts in the package can now be run from the command line with\n\n    oat <script_name>\n\nimporting functions in your own scripts can be done with\n\n    from oxDNA_analysis_tools.<file name> import <function>")

print("")

print("If you would like to enable autocomplete in Bash, either add\n\n    source oat-completion.sh\n\nto your .bashrc or run \n\n    sudo cp oat-completion.sh /etc/bash_completion.d/\n\nto make autocompletes available to all users.")
