import pandas as pd
import os
import argparse

# Create argument parser
parser = argparse.ArgumentParser(description="Split a CSV file into batches by group.")
parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to input CSV file")
parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory")
parser.add_argument("-b", "--batch_size", type=int, default=50, help="Number of groups per batch (default: 50)")

# Parse arguments
args = parser.parse_args()
input_file = args.input_file
output_dir = args.output_dir
batch_size = args.batch_size

# Read input CSV
data = pd.read_csv(input_file, sep=',', header=None)

# Group by the 5th column
grouped = data.groupby(data.iloc[:, 5])

# Get list of group names
group_names = list(grouped.groups.keys())

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Process in batches
for i in range(0, len(group_names), batch_size):
    current_batch = group_names[i:i + batch_size]
    current_data = pd.concat([grouped.get_group(name) for name in current_batch])

    output_file = os.path.join(output_dir, f"output_{i // batch_size + 1}.csv")
    current_data.to_csv(output_file, index=False, header=False)

    print(f"Created: {output_file}")
