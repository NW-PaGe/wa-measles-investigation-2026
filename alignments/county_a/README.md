# Multiple domestic measles importations within a two-week period, Washington State, January 2026

## Alignments for Washington Sequences and Their Matches - County A

The genome from County A ([PP_004J958](https://pathoplexus.org/seq/PP_004J958.2)) matched 161 sequences from South Carolina, so getting these matches into a post-hoc alignment was more involved (not practical to click through the tree in auspice to pull those data).

This assumes you have already run the phylogenetic workflow. If not, see the top-level README.md and run it, so that the dataset and tree are available.

To get the matching sequences, navigate to the `tree_fig` directory and activate the virtual environment:

```bash
source .env/bin/activate
```

Then navigate back to this directory and run `get_co_a_matches.py`. This pulls the matches with the County A genome and outputs the Pathoplexus accessions to a text file.

```bash
python3 get_co_a_matches.py
```

Filter the nextstrain dataset using the output from `get_co_a_matches.py`:

```
augur filter \
    --sequences ../../phylogenetic/results/sequences.fasta \
    --metadata ../../phylogenetic/results/metadata.tsv \
    --metadata-id-columns 'accession' \
    --exclude-all \
    --include county_a_matches.txt \
    --output-sequences county_a_filtered.fasta \
    --output-metadata county_a_filtered.tsv
```

Then align using `augur align`:

```bash
augur align --sequences county_a_filtered.fasta -o county_a_aligned.fasta
```

This alignment and it's statistics can be visualized by using the alignment viewer of your choice. We examined this alignment using a [web-based viewer](https://github.com/DOH-PNT0303/streamlit_align_and_view_app).
