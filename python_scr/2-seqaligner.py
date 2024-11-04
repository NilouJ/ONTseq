import os
import subprocess

# Define paths for reference, query directory, and output
ref_fasta_path = "/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye/5k/barcode03/assembly.fasta"
query_fastq_dir = "/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/ONTseq_Data/basecalling/pass/barcode03"
output_base = "/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_seqaligner"

# Ensure output directory exists
os.makedirs(output_base, exist_ok=True)

# Define output filenames based on the reference file
output_name = os.path.splitext(os.path.basename(ref_fasta_path))[0]
sam_output_path = os.path.join(output_base, f"{output_name}.sam")
bam_output_path = os.path.join(output_base, f"{output_name}.bam")

# Get all FASTQ files in the query directory
query_fastq_files = [os.path.join(query_fastq_dir, f) for f in os.listdir(query_fastq_dir) if f.endswith('.fastq.gz')]

# Step 1: Run minimap2 to align and produce a SAM file
print("Running minimap2 for alignment...")
minimap2_command = [
    "minimap2", "-ax", "map-ont", ref_fasta_path, *query_fastq_files, "-o", sam_output_path
]
try:
    subprocess.run(minimap2_command, check=True)
    print(f"SAM file generated at: {sam_output_path}")
except subprocess.CalledProcessError as e:
    print(f"Error running minimap2: {e}")
    exit(1)

# Step 2: Convert SAM to BAM and sort it
print("Converting SAM to BAM and sorting...")
samtools_view_command = ["samtools", "view", "-Sb", sam_output_path]
samtools_sort_command = ["samtools", "sort", "-o", bam_output_path]
try:
    # Pipe samtools view output directly into samtools sort
    view_process = subprocess.Popen(samtools_view_command, stdout=subprocess.PIPE)
    sort_process = subprocess.run(samtools_sort_command, stdin=view_process.stdout, check=True)
    view_process.stdout.close()
    print(f"BAM file generated at: {bam_output_path}")
except subprocess.CalledProcessError as e:
    print(f"Error converting SAM to BAM: {e}")
    exit(1)

# Step 3: Index the sorted BAM file
print("Indexing BAM file...")
samtools_index_command = ["samtools", "index", bam_output_path]
try:
    subprocess.run(samtools_index_command, check=True)
    print(f"Index file generated at: {bam_output_path}.bai")
except subprocess.CalledProcessError as e:
    print(f"Error indexing BAM file: {e}")
    exit(1)

print("Alignment, sorting, and indexing complete.")
