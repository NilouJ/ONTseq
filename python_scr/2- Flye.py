import subprocess
import os


def run_flye(fastq_dir, output_dir, genome_size, threads=6):
    # Verify directories exist
    if not os.path.isdir(fastq_dir):
        print(f"Error: The specified fastq directory '{fastq_dir}' does not exist.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory '{output_dir}' created.")

    # Construct the Flye command
    flye_cmd = [
        "flye",
        "--nano-raw", os.path.join(fastq_dir, "*.fastq.gz"),
        "--out-dir", output_dir,
        "--genome-size", genome_size,
        "--threads", str(threads)
    ]

    # Run the Flye command
    try:
        subprocess.run(flye_cmd, check=True)
        print("Flye assembly completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Flye: {e}")
    except FileNotFoundError:
        print("Error: Flye is not installed or not found in PATH.")


# Example usage:
# Customize these variables with your paths and parameters
fastq_directory = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/ONTseq_Data"
                   "/basecalling/pass/barcode01")
output_directory = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye/5k"
                    "/barcode01")
genome_size = "6.25k"

# Run the function
run_flye(fastq_directory, output_directory, genome_size)
