# Multiple domestic measles importations within a two-week period, Washington State, January 2026

## Data and Code for Reproducing the Analyses

### Dependencies:
These analyses can be run in Linux environments and Mac OS (untested).

Dependencies for Nextstrain build:
 -  Docker (https://docs.docker.com/engine/install/)
 -  Nextstrain CLI (https://docs.nextstrain.org/en/latest/install.html#install-nextstrain-cli)

Dependencies for phylogenetic tree manipulation/vizualization:
 
 -  uv (https://docs.astral.sh/uv/getting-started/installation/#installation-methods)

 ### Running the phylogenetic workflow
 
 These analyses were performed using the docker runtime for Nextstrain. Pull the same version:
 ```
 docker pull nextstrain/base:build-20260210T230050Z $$
 nextstrain setup --set-default docker
 ```

 Then navigate to the phylogenetic folder and run the workflow:

 ```
 cd phylogenetic &&
 nextstrain build --image nextstrain/base:build-20260210T230050Z . --configfile build-configs/wa_outbreak.yaml
 ```

 ### Process the tree with python

 ```
 cd tree_fig/ &&
 uv run main.py
 ```