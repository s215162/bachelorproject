import numpy as np
import pandas as pd
import os
import logging
import sys

def setup_logging(peptide):
    """Set up logging for the script."""
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/binder_processing_{peptide}.log"
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO,
    )
    
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)
    
    return logging

def process_cdr_files(peptide):
    """Process CDR files to add binder information and create filtered versions."""
    # Paths
    csv_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
    input_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    
    # Load CSV file
    logging.info(f"Loading CSV file: {csv_path}")
    df = pd.read_csv(csv_path)
    
    # Filter CSV for the specific peptide
    df_filtered = df[df['peptide_x'] == peptide]
    
    # Create binder lookup dictionary for faster access
    binder_dict = dict(zip(df_filtered['raw_index'], df_filtered['binder']))
    
    # Process each CDR file
    cdr_regions = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR1b', 'CDR2b', 'CDR3b']
    
    for cdr_region in cdr_regions:
        input_file = os.path.join(input_dir, f"{peptide}_{cdr_region}_embeddings.npz")
        
        if not os.path.exists(input_file):
            logging.warning(f"File not found: {input_file}")
            continue
            
        logging.info(f"Processing {cdr_region}")
        
        try:
            # Load the embeddings file
            data = np.load(input_file, allow_pickle=True)
            embeddings = data['embeddings']
            raw_indexes = data['raw_indexes']
            
            # Get binder values for each embedding
            binder_values = []
            filtered_embeddings = []
            filtered_raw_indexes = []
            
            for i, raw_index in enumerate(raw_indexes):
                binder_value = binder_dict.get(raw_index)
                if binder_value is not None:
                    binder_values.append(binder_value)
                    
                    # If this is a binder (value=1), add to filtered lists
                    if binder_value == 1:
                        filtered_embeddings.append(embeddings[i])
                        filtered_raw_indexes.append(raw_index)
                else:
                    logging.warning(f"No binder value found for raw_index {raw_index}")
            
            # Save full version with binder information
            full_output = os.path.join(input_dir, f"{peptide}_{cdr_region}_full.npz")
            np.savez(full_output, 
                     embeddings=embeddings,
                     raw_indexes=raw_indexes,
                     binder_values=binder_values)
            logging.info(f"Saved full data with binder values to {full_output}")
            
            # Save binder-only version
            if filtered_embeddings:
                binder_output = os.path.join(input_dir, f"{peptide}_{cdr_region}_binder.npz")
                np.savez(binder_output,
                         embeddings=np.array(filtered_embeddings, dtype=object),
                         raw_indexes=filtered_raw_indexes)
                logging.info(f"Saved binder-only data to {binder_output}")
                logging.info(f"Number of binders: {len(filtered_embeddings)}")
            else:
                logging.warning(f"No binders found for {cdr_region}")
                
        except Exception as e:
            logging.error(f"Error processing {cdr_region}: {e}")
            continue

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptide = sys.argv[1]
    
    # Setup logging
    logger = setup_logging(peptide)
    
    try:
        process_cdr_files(peptide)
        logging.info(f"Processing completed for peptide {peptide}")
    
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()