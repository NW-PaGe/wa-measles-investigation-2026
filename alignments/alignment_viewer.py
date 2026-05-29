#!/usr/bin/env python3
"""
Browser-based sequence alignment viewer with alignment capability.
"""

import streamlit as st
import tempfile
import os
import subprocess
from pathlib import Path
from collections import defaultdict

def read_aligned_fasta_from_string(fasta_content):
    """Read aligned FASTA from string."""
    sequences = {}
    current_name = None
    current_seq = []

    for line in fasta_content.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            if current_name:
                sequences[current_name] = ''.join(current_seq)
            current_name = line[1:]
            current_seq = []
        else:
            current_seq.append(line)

    if current_name:
        sequences[current_name] = ''.join(current_seq)

    return sequences

def check_if_aligned(sequences):
    """Check if sequences are already aligned (all same length)."""
    if not sequences:
        return False
    lengths = [len(seq) for seq in sequences.values()]
    return len(set(lengths)) == 1

def align_sequences_mafft(fasta_content):
    """Align sequences using MAFFT."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        input_file = f.name

    output_file = input_file + '.aligned'

    try:
        # Run MAFFT
        result = subprocess.run(
            ['mafft', '--auto', '--thread', '-1', input_file],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        if result.returncode != 0:
            raise Exception(f"MAFFT failed: {result.stderr}")

        aligned_content = result.stdout

        # Cleanup
        os.unlink(input_file)

        return aligned_content, None

    except FileNotFoundError:
        return None, "MAFFT not found. Please install MAFFT first."
    except subprocess.TimeoutExpired:
        return None, "Alignment timed out (>5 minutes). Try with fewer/shorter sequences."
    except Exception as e:
        return None, str(e)
    finally:
        if os.path.exists(input_file):
            os.unlink(input_file)

def align_sequences_muscle(fasta_content):
    """Align sequences using MUSCLE."""
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(fasta_content)
        input_file = f.name

    output_file = input_file + '.aligned'

    try:
        # Run MUSCLE (v5 syntax)
        result = subprocess.run(
            ['muscle', '-align', input_file, '-output', output_file],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            # Try MUSCLE v3 syntax
            result = subprocess.run(
                ['muscle', '-in', input_file, '-out', output_file],
                capture_output=True,
                text=True,
                timeout=300
            )

            if result.returncode != 0:
                raise Exception(f"MUSCLE failed: {result.stderr}")

        with open(output_file, 'r') as f:
            aligned_content = f.read()

        # Cleanup
        os.unlink(input_file)
        os.unlink(output_file)

        return aligned_content, None

    except FileNotFoundError:
        return None, "MUSCLE not found. Please install MUSCLE first."
    except subprocess.TimeoutExpired:
        return None, "Alignment timed out (>5 minutes). Try with fewer/shorter sequences."
    except Exception as e:
        return None, str(e)
    finally:
        if os.path.exists(input_file):
            os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

def calculate_pairwise_distances(sequences):
    """Calculate pairwise SNP distances between all sequence pairs (non-redundant)."""
    seq_names = list(sequences.keys())
    ambiguity_codes = set('RYSWKMBDHVN')

    pairwise_distances = []
    distance_matrix = {}  # For heatmap

    # Initialize matrix
    for name in seq_names:
        distance_matrix[name] = {}

    for i in range(len(seq_names)):
        for j in range(i + 1, len(seq_names)):
            name1 = seq_names[i]
            name2 = seq_names[j]
            seq1 = sequences[name1]
            seq2 = sequences[name2]

            snps = 0
            usable_positions = 0

            for pos in range(len(seq1)):
                base1 = seq1[pos].upper()
                base2 = seq2[pos].upper()

                # Skip if either has a gap
                if base1 == '-' or base2 == '-':
                    continue

                # Skip if either has ambiguity code
                if base1 in ambiguity_codes or base2 in ambiguity_codes:
                    continue

                usable_positions += 1

                if base1 != base2:
                    snps += 1

            identity = 100 * (1 - snps / usable_positions) if usable_positions > 0 else 0

            pairwise_distances.append({
                'seq1': name1,
                'seq2': name2,
                'snps': snps,
                'usable_positions': usable_positions,
                'identity': identity
            })

            # Store in matrix (both directions for easy lookup)
            distance_matrix[name1][name2] = snps
            distance_matrix[name2][name1] = snps

    return pairwise_distances, distance_matrix, seq_names

def find_differences(sequences):
    """Find positions where sequences differ (ignoring gaps and ambiguity codes)."""
    if not sequences:
        return [], [], []

    seq_list = list(sequences.values())
    align_len = len(seq_list[0])

    ambiguity_codes = set('RYSWKMBDHVN')

    differences = []
    gap_positions = []
    ambiguous_positions = []

    for pos in range(align_len):
        bases = [seq[pos].upper() for seq in seq_list]
        bases_set = set(bases)

        has_gap = '-' in bases_set
        has_ambiguity = bool(bases_set & ambiguity_codes)

        if has_gap:
            gap_positions.append(pos)
        elif has_ambiguity:
            ambiguous_positions.append(pos)
        else:
            if len(bases_set) > 1:
                differences.append(pos)

    return differences, gap_positions, ambiguous_positions

def generate_html_visualization(sequences, differences, gap_positions, ambiguous_positions, pairwise_distances, distance_matrix, seq_names):
    """Generate HTML visualization as string."""
    align_len = len(list(sequences.values())[0])

    # Calculate min and max SNPs for color scaling
    all_snps = [pair['snps'] for pair in pairwise_distances] if pairwise_distances else [0]
    min_snps = min(all_snps) if all_snps else 0
    max_snps = max(all_snps) if all_snps else 0

    html = """
<!DOCTYPE html>
<html>
<head>
    <title>Alignment Viewer</title>
    <style>
        body {
            font-family: 'Courier New', monospace;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
        }
        .stats {
            background-color: white;
            padding: 15px;
            margin-bottom: 20px;
            border-radius: 5px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .heatmap-container {
            overflow-x: auto;
            margin-top: 15px;
        }
        .heatmap-table {
            border-collapse: collapse;
            margin: 0 auto;
        }
        .heatmap-table th,
        .heatmap-table td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: center;
            min-width: 80px;
            font-size: 11px;
        }
        .heatmap-table th {
            background-color: #34495e;
            color: white;
            font-weight: bold;
            writing-mode: vertical-rl;
            text-orientation: mixed;
            max-width: 30px;
            padding: 8px 4px;
        }
        .heatmap-table td.row-header {
            background-color: #34495e;
            color: white;
            font-weight: bold;
            text-align: left;
            writing-mode: horizontal-tb;
            max-width: 200px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }
        .heatmap-table td.heatmap-cell {
            font-weight: bold;
            cursor: help;
        }
        .heatmap-table td.diagonal {
            background-color: #95a5a6;
            color: #7f8c8d;
        }
        .pairwise-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }
        .pairwise-table th,
        .pairwise-table td {
            border: 1px solid #ddd;
            padding: 8px;
            text-align: left;
        }
        .pairwise-table th {
            background-color: #34495e;
            color: white;
            font-weight: bold;
        }
        .pairwise-table tr:nth-child(even) {
            background-color: #f2f2f2;
        }
        .pairwise-table tr:hover {
            background-color: #e8e8e8;
        }
        .alignment {
            background-color: white;
            padding: 15px;
            border-radius: 5px;
            overflow-x: auto;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .alignment-content {
            min-width: fit-content;
        }
        /* Styling for scrollbar */
        .alignment::-webkit-scrollbar {
            height: 12px;
        }
        .alignment::-webkit-scrollbar-track {
            background: #f1f1f1;
            border-radius: 10px;
        }
        .alignment::-webkit-scrollbar-thumb {
            background: #888;
            border-radius: 10px;
        }
        .alignment::-webkit-scrollbar-thumb:hover {
            background: #555;
        }
        .seq-name {
            display: inline-block;
            width: 300px;
            font-weight: bold;
            background-color: #ecf0f1;
            padding: 2px 5px;
            margin: 2px 0;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }
        .seq-line {
            margin: 2px 0;
            white-space: nowrap;
            display: flex;
            align-items: center;
        }
        .seq-sequence {
            flex-shrink: 0;
            white-space: nowrap;
        }
        .match {
            color: #27ae60;
        }
        .diff {
            background-color: #e74c3c;
            color: white;
            font-weight: bold;
            padding: 0 2px;
        }
        .gap {
            color: #95a5a6;
        }
        .ambig {
            background-color: #f39c12;
            color: white;
            font-weight: bold;
            padding: 0 2px;
        }
        .position-marker {
            color: #7f8c8d;
            font-size: 10px;
            margin-top: 15px;
        }
        .legend {
            color: #7f8c8d;
            margin-bottom: 15px;
        }
        .snp-range {
            font-size: 12px;
            color: #7f8c8d;
            margin-top: 5px;
        }
    </style>
</head>
<body>
    <div class="stats">
        <h2>Summary Statistics</h2>
"""

    html += f"        <p><strong>Number of sequences:</strong> {len(sequences)}</p>\n"
    html += f"        <p><strong>Alignment length:</strong> {align_len:,} bp</p>\n"
    html += f"        <p><strong>Positions with gaps:</strong> {len(gap_positions):,}</p>\n"
    html += f"        <p><strong>Positions with ambiguity codes:</strong> {len(ambiguous_positions):,}</p>\n"
    usable_positions = align_len - len(gap_positions) - len(ambiguous_positions)
    html += f"        <p><strong>Usable positions:</strong> {usable_positions:,}</p>\n"
    html += f"        <p><strong>Variable sites (any variation):</strong> {len(differences):,}</p>\n"

    # Add SNP range for multiple sequences
    if len(sequences) > 2 and pairwise_distances:
        html += f"        <p><strong>Pairwise SNP range:</strong> {min_snps:,} - {max_snps:,}</p>\n"
    elif len(sequences) == 2 and pairwise_distances:
        html += f"        <p><strong>Pairwise SNPs:</strong> {pairwise_distances[0]['snps']:,}</p>\n"
        html += f"        <p><strong>Identity:</strong> {pairwise_distances[0]['identity']:.2f}%</p>\n"

    # Add heatmap distance matrix if there are multiple sequences
    if len(sequences) > 1 and distance_matrix:
        html += """
        <h3>Pairwise SNP Distance Matrix (Heatmap)</h3>
        <p class="snp-range">Color scale: Blue (fewer SNPs) → Red (more SNPs)</p>
        <div class="heatmap-container">
            <table class="heatmap-table">
                <thead>
                    <tr>
                        <th></th>
"""
        # Column headers
        for name in seq_names:
            display_name = name[:20] + '...' if len(name) > 20 else name
            html += f"                        <th title='{name}'>{display_name}</th>\n"

        html += """                    </tr>
                </thead>
                <tbody>
"""

        # Data rows
        for i, name1 in enumerate(seq_names):
            display_name1 = name1[:30] + '...' if len(name1) > 30 else name1
            html += f"                    <tr>\n"
            html += f"                        <td class='row-header' title='{name1}'>{display_name1}</td>\n"

            for j, name2 in enumerate(seq_names):
                if i == j:
                    # Diagonal
                    html += f"                        <td class='diagonal heatmap-cell'>-</td>\n"
                else:
                    snp_count = distance_matrix[name1][name2]

                    # Calculate color based on SNP count (blue to red scale)
                    if max_snps > min_snps:
                        ratio = (snp_count - min_snps) / (max_snps - min_snps)
                    else:
                        ratio = 0

                    # Blue (low SNPs) to Red (high SNPs) via white
                    if ratio < 0.5:
                        # Blue to white
                        r = int(240 + (255 - 240) * (ratio * 2))
                        g = int(248 + (255 - 248) * (ratio * 2))
                        b = 255
                    else:
                        # White to red
                        r = 255
                        g = int(255 - (255 - 180) * ((ratio - 0.5) * 2))
                        b = int(255 - (255 - 180) * ((ratio - 0.5) * 2))

                    color = f"rgb({r},{g},{b})"
                    text_color = "#000" if ratio < 0.7 else "#fff"

                    html += f"                        <td class='heatmap-cell' style='background-color: {color}; color: {text_color}' title='{name1} vs {name2}: {snp_count:,} SNPs'>{snp_count:,}</td>\n"

            html += "                    </tr>\n"

        html += """                </tbody>
            </table>
        </div>
"""

        # Add detailed table
        html += """
        <h3>Detailed Pairwise Comparisons</h3>
        <table class="pairwise-table">
            <thead>
                <tr>
                    <th>Sequence 1</th>
                    <th>Sequence 2</th>
                    <th>SNPs</th>
                    <th>Usable Positions</th>
                    <th>Identity (%)</th>
                </tr>
            </thead>
            <tbody>
"""
        for pair in pairwise_distances:
            html += f"""                <tr>
                    <td>{pair['seq1']}</td>
                    <td>{pair['seq2']}</td>
                    <td>{pair['snps']:,}</td>
                    <td>{pair['usable_positions']:,}</td>
                    <td>{pair['identity']:.2f}%</td>
                </tr>
"""
        html += """            </tbody>
        </table>
"""

    html += """
    </div>

    <div class="alignment">
        <h2>Alignment Visualization</h2>
        <div class="legend">
            Legend: <span class="match">Match</span> |
            <span class="diff">SNP</span> |
            <span class="gap">Gap</span> |
            <span class="ambig">Ambiguous</span>
        </div>
        <div class="alignment-content">
"""

    diff_set = set(differences)
    gap_set = set(gap_positions)
    ambig_set = set(ambiguous_positions)

    chunk_size = 100
    for chunk_start in range(0, align_len, chunk_size):
        chunk_end = min(chunk_start + chunk_size, align_len)

        html += f"        <div style='margin: 20px 0;'>\n"
        html += f"        <div class='position-marker'>Position {chunk_start + 1}-{chunk_end}</div>\n"

        for name in seq_names:
            seq = sequences[name]
            chunk = seq[chunk_start:chunk_end]

            html += f"        <div class='seq-line'>\n"
            html += f"            <span class='seq-name' title='{name}'>{name[:50]}</span>\n"
            html += f"            <span class='seq-sequence'>"

            for i, base in enumerate(chunk):
                pos = chunk_start + i
                if pos in ambig_set:
                    html += f"<span class='ambig'>{base}</span>"
                elif pos in diff_set:
                    html += f"<span class='diff'>{base}</span>"
                elif base == '-':
                    html += f"<span class='gap'>{base}</span>"
                else:
                    html += f"<span class='match'>{base}</span>"

            html += "</span>\n        </div>\n"

        html += "        </div>\n"

    html += """
        </div>
    </div>
</body>
</html>
"""

    return html

# Streamlit App
st.set_page_config(page_title="Sequence Alignment Viewer", layout="wide")

st.title("Sequence Alignment Viewer")
st.markdown("Upload sequences to align and visualize differences")

# Sidebar options
with st.sidebar:
    st.header("Settings")

    alignment_tool = st.selectbox(
        "Alignment Tool",
        ["MAFFT (recommended)", "MUSCLE"],
        help="Choose alignment algorithm. MAFFT is faster and recommended for most cases."
    )

    auto_align = st.checkbox(
        "Auto-align if needed",
        value=True,
        help="Automatically align sequences if they're not already aligned"
    )

    st.markdown("---")
    st.markdown("""
    **Installation:**
    ```bash
    # MAFFT (recommended)
    conda install -c bioconda mafft
    # or
    sudo apt install mafft

    # MUSCLE (alternative)
    conda install -c bioconda muscle
    ```
    """)

# File uploader
uploaded_file = st.file_uploader(
    "Choose a FASTA file",
    type=['fasta', 'fa', 'fna', 'txt'],
    help="Upload aligned or unaligned FASTA sequences"
)

if uploaded_file is not None:
    # Read the file content
    content = uploaded_file.read().decode('utf-8')

    # Parse sequences
    try:
        sequences = read_aligned_fasta_from_string(content)

        if not sequences:
            st.error("No sequences found in file")
        elif len(sequences) < 2:
            st.error("Need at least 2 sequences to compare")
        else:
            # Check if sequences are aligned
            is_aligned = check_if_aligned(sequences)

            if not is_aligned:
                if auto_align:
                    st.info("Sequences are not aligned. Aligning now...")

                    # Choose alignment tool
                    if "MAFFT" in alignment_tool:
                        aligned_content, error = align_sequences_mafft(content)
                    else:
                        aligned_content, error = align_sequences_muscle(content)

                    if error:
                        st.error(f"Alignment failed: {error}")
                        st.stop()
                    else:
                        st.success("Alignment complete!")
                        sequences = read_aligned_fasta_from_string(aligned_content)

                        # Offer download of aligned sequences
                        st.download_button(
                            label="Download Aligned FASTA",
                            data=aligned_content,
                            file_name="aligned_sequences.fasta",
                            mime="text/plain"
                        )
                else:
                    st.warning("Sequences are not aligned. Enable 'Auto-align if needed' in settings to align them.")
                    st.stop()
            else:
                st.success("Sequences are already aligned")

            # Find differences
            differences, gap_positions, ambiguous_positions = find_differences(sequences)

            # Calculate pairwise distances
            pairwise_distances, distance_matrix, seq_names_ordered = calculate_pairwise_distances(sequences)

            # Display summary in sidebar
            with st.sidebar:
                st.header("Quick Stats")
                align_len = len(list(sequences.values())[0])
                usable_positions = align_len - len(gap_positions) - len(ambiguous_positions)

                st.metric("Sequences", len(sequences))
                st.metric("Alignment Length", f"{align_len:,} bp")

                if len(sequences) == 2 and pairwise_distances:
                    st.metric("SNPs", pairwise_distances[0]['snps'])
                    st.metric("Identity", f"{pairwise_distances[0]['identity']:.2f}%")
                elif len(sequences) > 2 and pairwise_distances:
                    min_snps = min(p['snps'] for p in pairwise_distances)
                    max_snps = max(p['snps'] for p in pairwise_distances)
                    st.metric("SNP Range", f"{min_snps} - {max_snps}")

            # Generate and display HTML
            html_content = generate_html_visualization(
                sequences, differences, gap_positions, ambiguous_positions,
                pairwise_distances, distance_matrix, seq_names_ordered
            )

            # Display in iframe
            st.components.v1.html(html_content, height=800, scrolling=True)

            # Download button
            st.download_button(
                label="Download HTML Visualization",
                data=html_content,
                file_name="alignment_visualization.html",
                mime="text/html"
            )

    except Exception as e:
        st.error(f"Error processing file: {str(e)}")
        st.exception(e)

else:
    # Show example
    st.info("Upload a FASTA file to get started")

    with st.expander("About this tool"):
        st.markdown("""
        This tool can both **align** and **visualize** sequences:

        **Features:**
        - Auto-align unaligned sequences using MAFFT or MUSCLE
        - Visualize differences with color coding
        - Calculate identity and SNP statistics
        - Download aligned sequences and HTML visualizations

        **Visualization:**
        - **SNPs** are highlighted in red
        - **Gaps** are shown in gray
        - **Ambiguous bases** (IUPAC codes) are shown in orange
        - **Matches** are shown in green

        **SNP Counting:**
        - Pairwise distances calculated for all sequence pairs (non-redundant)
        - Matches `snp-dists` behavior for pairwise comparisons
        - Excludes positions with gaps in either sequence
        - Excludes positions with ambiguity codes (RYSWKMBDHVN)
        - Only counts true nucleotide differences
        - For multiple sequences, displays distance matrix table

        **Supported formats:**
        - Aligned or unaligned FASTA files
        - Multiple sequences for comparison
        """)

    with st.expander("Example usage"):
        st.markdown("""
        1. Upload your FASTA file (aligned or unaligned)
        2. If unaligned, sequences will be aligned automatically
        3. View the interactive visualization
        4. Download the results

        Works great for:
        - Comparing viral genomes
        - Checking sequence variants
        - Quality control of sequencing data
        - Identifying SNPs and mutations
        """)
