# transpose_tsv_native.py
import sys

def transpose_tsv(input_file, output_file):
    """
    Reads a tab-separated matrix from input_file,
    transposes it natively, and writes it to output_file.
    """
    try:
        # 1. Import: Read and parse the tab-separated file
        with open(input_file, 'r', encoding='utf-8') as f:
            # Split lines by tab, strip trailing newlines/spaces
            matrix = [line.strip().split('\t') for line in f if line.strip()]

        if not matrix:
            print(f"Error: {input_file} is empty.")
            return

        # 2. Transpose: Use native zip unpacking
        transposed = list(zip(*matrix))

        # 3. Export: Write the transposed matrix back as tab-separated
        with open(output_file, 'w', encoding='utf-8') as f:
            for row in transposed:
                f.write('\t'.join(row) + '\n')

        print(f"Successfully transposed {input_file} -> {output_file}")

    except FileNotFoundError:
        print(f"Error: The file {input_file} was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    # Set default filenames (can be modified as needed)
    INPUT_FILENAME = 'weighted.inp'
    OUTPUT_FILENAME = 'fuzzy-rank-energy.inp'

    transpose_tsv(INPUT_FILENAME, OUTPUT_FILENAME)
