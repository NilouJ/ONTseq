# import pysam
# import plotly.graph_objects as go
# from plotly.subplots import make_subplots
#
# # Path to the BAM file with MD tag
# bam_path = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/results_minimap"
#             "/barcode03/aligned_reads_filtered_md.bam")
# ref_fasta_path = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye/50k"
#                   "/barcode03/assembly.fasta")
# samfile = pysam.AlignmentFile(bam_path, "rb")
# fasta = pysam.FastaFile(ref_fasta_path)
#
# # Initial parameters for the window
# window_size = 200
# ref_seq = fasta.fetch("contig_1")
# total_length = len(ref_seq)
#
# # Define color scheme for nucleotides
# colors = {'A': '#1f77b4', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'}
#
#
# def get_trace_data(start_pos, window_size):
#     end_pos = min(start_pos + window_size, total_length)
#     x_values = list(range(start_pos, end_pos))
#
#     # Initialize containers for nucleotides
#     nucleotide_confidence = {'A': [], 'T': [], 'C': [], 'G': []}
#     nucleotide_positions = {'A': [], 'T': [], 'C': [], 'G': []}
#
#     # Extract alignment info for this window
#     for read in samfile.fetch("contig_1", start=start_pos, end=end_pos):
#         for query_pos, ref_pos, base in read.get_aligned_pairs(matches_only=True, with_seq=True):
#             if ref_pos is not None and start_pos <= ref_pos < end_pos and base in nucleotide_confidence:
#                 confidence = read.query_qualities[query_pos] if read.query_qualities else 30
#                 nucleotide_confidence[base].append(confidence)
#                 nucleotide_positions[base].append(ref_pos)
#
#     # Create traces
#     traces = [
#         go.Bar(
#             x=nucleotide_positions[nucleotide],
#             y=nucleotide_confidence[nucleotide],
#             name=f"{nucleotide} Confidence",
#             marker_color=colors[nucleotide],
#             showlegend=True
#         ) for nucleotide in ['A', 'T', 'C', 'G']
#     ]
#
#     # Reference sequence labels
#     ref_labels = go.Scatter(
#         x=x_values,
#         y=[1] * len(x_values),
#         text=[ref_seq[i] for i in x_values],
#         mode="text",
#         textfont=dict(size=10, color='gray'),
#         showlegend=False
#     )
#
#     return traces, ref_labels
#
#
# # Setup initial figure with first window
# initial_start = 0
# traces, ref_labels = get_trace_data(initial_start, window_size)
#
# # Initialize figure
# fig = make_subplots(
#     rows=2, cols=1, shared_xaxes=True,
#     row_heights=[0.85, 0.15],
#     vertical_spacing=0.03
# )
#
# # Add initial traces
# for trace in traces:
#     fig.add_trace(trace, row=1, col=1)
# fig.add_trace(ref_labels, row=2, col=1)
#
# # Configure layout with slider
# fig.update_layout(
#     title="Nucleotide Confidence Levels (Sliding Window of 200)",
#     xaxis_title="Position",
#     yaxis_title="Confidence Level",
#     height=500,
#     template="plotly_white",
#     font=dict(family="Arial", size=12),
#     sliders=[{
#         "active": 0,
#         "currentvalue": {"prefix": "Position: ", "font": {"size": 12}},
#         "pad": {"t": 50},
#         "steps": [{
#             "label": f"{start}-{start + window_size}",
#             "method": "relayout",
#             "args": [{"xaxis.range": [start, start + window_size]}]
#         } for start in range(0, total_length, window_size)]
#     }]
# )
#
# # Update style and show plot
# fig.update_xaxes(matches='x', showline=True, linewidth=1, linecolor='black')
# fig.update_yaxes(title_text="Confidence Level", showline=True, linewidth=1, linecolor='black', row=1, col=1)
# fig.update_yaxes(showticklabels=False, row=2, col=1)
# fig.show()

import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import plotly.graph_objs as go
import pysam

# Path to the BAM file with MD tag
bam_path = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/results_minimap"
            "/barcode03/aligned_reads_filtered_md.bam")
ref_fasta_path = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye/50k"
                  "/barcode03/assembly.fasta")

# Load BAM and reference FASTA files
samfile = pysam.AlignmentFile(bam_path, "rb")
fasta = pysam.FastaFile(ref_fasta_path)

# Initialize Dash app
app = dash.Dash(__name__)

# Configure layout
app.layout = html.Div([
    html.H1("ONT Nucleotide Confidence", style={"text-align": "center"}),
    dcc.Graph(id="nucleotide-graph"),
    dcc.RangeSlider(
        id="position-slider",
        min=0,
        max=len(fasta.fetch("contig_1")) - 200,
        value=[0, 200],
        step=100,
        marks={i: str(i) for i in range(0, len(fasta.fetch("contig_1")), 1000)},
        tooltip={"placement": "bottom", "always_visible": True}
    ),
    html.Div(id="slider-output-container", style={"text-align": "center"})
])

# Define color scheme for nucleotides
colors = {'A': '#1f77b4', 'T': '#ff7f0e', 'C': '#2ca02c', 'G': '#d62728'}


# Callback to update graph based on slider position
@app.callback(
    Output("nucleotide-graph", "figure"),
    [Input("position-slider", "value")]
)
def update_graph(range_vals):
    start_pos, end_pos = range_vals
    x_values = list(range(start_pos, end_pos))

    # Initialize containers for nucleotides
    nucleotide_confidence = {"A": [], "T": [], "C": [], "G": []}
    nucleotide_positions = {"A": [], "T": [], "C": [], "G": []}

    # Extract alignment info for this window
    for read in samfile.fetch("contig_1", start=start_pos, end=end_pos):
        for query_pos, ref_pos, base in read.get_aligned_pairs(matches_only=True, with_seq=True):
            if ref_pos is not None and start_pos <= ref_pos < end_pos and base in nucleotide_confidence:
                confidence = read.query_qualities[query_pos] if read.query_qualities else 30
                nucleotide_confidence[base].append(confidence)
                nucleotide_positions[base].append(ref_pos)

    # Create traces for each nucleotide
    traces = [
        go.Scatter(
            x=nucleotide_positions[nucleotide],
            y=nucleotide_confidence[nucleotide],
            mode="markers",
            marker=dict(color=colors[nucleotide], size=6),
            name=f"{nucleotide} Confidence"
        ) for nucleotide in ["A", "T", "C", "G"]
    ]

    # Construct figure
    figure = {
        "data": traces,
        "layout": go.Layout(
            title=f"Nucleotide Confidence Levels for Positions {start_pos}-{end_pos}",
            xaxis={"title": "Position", "range": [start_pos, end_pos]},
            yaxis={"title": "Confidence Level"},
            showlegend=True,
            margin=dict(l=40, r=40, t=40, b=40),
            template="plotly_white",
            hovermode="closest"
        )
    }

    return figure


# Run the app
if __name__ == "__main__":
    app.run_server(debug=True)
