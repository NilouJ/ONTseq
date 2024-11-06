import subprocess
import os
import glob


def run_flye(fastq_dir, output_dir, genome_size, threads=6, use_meta=False):
    if not os.path.isdir(fastq_dir):
        print(f"Error: The specified fastq directory '{fastq_dir}' does not exist.")
        return

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Output directory '{output_dir}' created.")

    fastq_files = glob.glob(os.path.join(fastq_dir, "*.fastq.gz"))
    if not fastq_files:
        print(f"Error: No .fastq.gz files found in the specified directory '{fastq_dir}'.")
        return

    flye_cmd = [
                   "flye",
                   "--nano-raw"
               ] + fastq_files + [
                   "--out-dir", output_dir,
                   "--genome-size", genome_size,
                   "--threads", str(threads)
               ]

    if use_meta:
        flye_cmd.append("--meta")

    # Print the constructed command for comparison
    print("Generated command:", " ".join(flye_cmd))

    try:
        subprocess.run(flye_cmd, check=True)
        print("Flye assembly completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running Flye: {e}")
    except FileNotFoundError:
        print("Error: Flye is not installed or not found in PATH.")


# Example usage
fastq_directory = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/ONTseq_Data"
                   "/basecalling/pass/barcode03")
output_directory = ("/Users/MQ10005295/Library/CloudStorage/OneDrive-MacquarieUniversity/Nilou/Nanopore/result_flye/6k"
                    "/barcode03")
genome_size = "6221"

run_flye(fastq_directory, output_directory, genome_size, threads=6, use_meta=True)
