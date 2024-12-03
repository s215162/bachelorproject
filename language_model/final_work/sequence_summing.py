import numpy as np
import pandas as pd
import os
import logging
import sys

def setup_logging(peptide):
    """Set up logging for the script."""
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/embedding_extraction_{peptide}.log"
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

def find_cdr_positions(tcr_full, cdr_seq):
    """Find the start and end positions of a CDR sequence in the full TCR sequence."""
    start_pos = tcr_full.find(cdr_seq)
    if start_pos == -1:
        return None, None
    end_pos = start_pos + len(cdr_seq)
    return start_pos, end_pos

def extract_embeddings(peptide):
    """Extract embeddings based on CDR regions."""
    # Paths
    npz_path = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/final_esm_embedding_matrix_full_data_{peptide}.npz"
    csv_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
    output_dir = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/cdr_embeddings_{peptide}"
    
    os.makedirs(output_dir, exist_ok=True)
    
    logging.info(f"Loading NPZ file: {npz_path}")
    npz_data = np.load(npz_path, allow_pickle=True)
    embeddings = npz_data['embeddings']
    raw_indexes = npz_data['raw_indexes']
    
    logging.info(f"Loading CSV file: {csv_path}")
    df = pd.read_csv(csv_path)
    
    cdr_regions = ['CDR1a', 'CDR2a', 'CDR3a', 'CDR1b', 'CDR2b', 'CDR3b']
    
    for cdr_region in cdr_regions:
        logging.info(f"Processing {cdr_region}")
        
        df_filtered = df[
            (df['peptide_x'] == peptide) & 
            (df[cdr_region].notna()) & 
            (df[cdr_region] != '')
        ]
        
        if len(df_filtered) == 0:
            logging.warning(f"No rows found for {cdr_region} in peptide {peptide}")
            continue
        
        region_embeddings = []
        region_raw_indexes = []
        
        for _, row in df_filtered.iterrows():
            try:
                idx = raw_indexes.tolist().index(row['raw_index'])
                full_embedding = embeddings[idx]
                
                # Find the CDR sequence position in the full TCR sequence
                start_pos, end_pos = find_cdr_positions(row['tcr_full'], row[cdr_region])
                
                if start_pos is not None:
                    # Extract only the rows corresponding to the CDR sequence
                    # Add 1 to account for the start token in ESM embeddings
                    cdr_embedding = full_embedding[start_pos+1:end_pos+1]
                    
                    logging.info(f"CDR sequence length: {len(row[cdr_region])}, Embedding shape: {cdr_embedding.shape}")
                    
                    region_embeddings.append(cdr_embedding)
                    region_raw_indexes.append(row['raw_index'])
                else:
                    logging.warning(f"Could not find CDR sequence in full TCR for raw_index {row['raw_index']}")
                    
            except ValueError:
                logging.warning(f"Raw index {row['raw_index']} not found in embeddings")
        
        if region_embeddings:
            region_embeddings_array = np.array(region_embeddings, dtype=object)
            
            output_file = os.path.join(output_dir, f"{peptide}_{cdr_region}_embeddings.npz")
            np.savez(output_file, 
                     embeddings=region_embeddings_array, 
                     raw_indexes=region_raw_indexes)
            
            logging.info(f"Saved embeddings for {cdr_region} to {output_file}")
            logging.info(f"Number of embeddings: {len(region_embeddings)}")
        else:
            logging.warning(f"No embeddings found for {cdr_region}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <peptide>")
        sys.exit(1)
    
    peptide = sys.argv[1]
    setup_logging(peptide)
    
    try:
        extract_embeddings(peptide)
        logging.info(f"Processing completed for peptide {peptide}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()