import pandas as pd

peptides = ["ELAGIGILTV", "GILGFVFTL", "GLCTLVAML", "LLWNGPMAV", "RAKFKQLL", "YLQPRTFLL"]
weights = ["CDR1a", "CDR1b", "CDR2a", "CDR2b", "CDR3a", "CDR3b", "TCRa", "TCRb"]
set = ["binders", "swaps"]

base_directory = "/net/mimer/mnt/tank/projects2/emison/language_model/27th_Oct_new_matrices/04.11_new_esm/12.11/"

for peptide in peptides:
    for weight in weights:
        for _type in set:
            file_path = base_directory + f"sequence_summed_vectors_{_type}_{peptide}_{weight}.csv"

            try:
                # Load the dataset into a pandas DataFrame
                df = pd.read_csv(file_path)

                # Function to remove 'tensor(' and ')' and convert to float
                def convert_tensor_to_float(value):
                    if isinstance(value, str) and value.startswith('tensor(') and value.endswith(')'):
                        return float(value[7:-1])  # Extract number from 'tensor(...)'
                    return value  # If it's not a tensor string, return the original value

                # Apply the conversion function to each column except the first one (index)
                for column in df.columns[1:]:
                    df[column] = df[column].map(convert_tensor_to_float)

                # Save the modified dataset to a new CSV file
                output_path = base_directory + f"{weight.lower()}/sequence_summed_vectors_{_type}_{peptide}_{weight}_np.csv"
                df.to_csv(output_path, index=False)

                print(f"Dataset has been processed and saved to {output_path}")

            except Exception as e:
                print(f"An error occurred while processing {file_path}: {e}")

