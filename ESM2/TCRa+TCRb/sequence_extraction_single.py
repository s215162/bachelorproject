import numpy as np
import pandas as pd
import os
import logging
import sys

def setup_logging(peptide):
   """Set up logging for the script."""
   # Update paths to fit
   log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/embedding_extraction_single_chain_{peptide}.log"
   # logging configuration
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

def find_cdr_positions(tcr_chain, cdr_seq):
   """Find the start and end positions of a CDR sequence in the TCR chain sequence.
   
   Returns None, None if CDR sequence not found in chain.
   """
   # Search for CDR sequence within the full chain
   start_pos = tcr_chain.find(cdr_seq)
   if start_pos == -1:
       return None, None
   end_pos = start_pos + len(cdr_seq)
   return start_pos, end_pos

def extract_embeddings_for_chain(peptide, chain_type, cdr_regions):
   """Extract embeddings for specific CDR regions from a TCR chain.
   
   Args:
       peptide: Target peptide identifier
       chain_type: Either 'TCRa' or 'TCRb'
       cdr_regions: List of CDR regions to process (e.g. ['CDR1a', 'CDR2a', 'CDR3a'])
   """
   # Define file paths for embeddings and sequence data. Update based on user directory
   npz_path = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/final_esm_embedding_matrix_{chain_type}_data_{peptide}.npz"
   csv_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
   output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
   
   os.makedirs(output_dir, exist_ok=True)
   
   # Load pre-computed ESM embeddings and sequence data
   logging.info(f"Loading NPZ file for {chain_type}: {npz_path}")
   npz_data = np.load(npz_path, allow_pickle=True)
   embeddings = npz_data['embeddings']
   raw_indexes = npz_data['raw_indexes']
   
   logging.info(f"Loading CSV file: {csv_path}")
   df = pd.read_csv(csv_path)
   
   # Process each CDR region separately
   for cdr_region in cdr_regions:
       logging.info(f"Processing {cdr_region}")
       
       # Filter sequences for current peptide and CDR region
       df_filtered = df[
           (df['peptide_x'] == peptide) & 
           (df[cdr_region].notna()) & 
           (df[cdr_region] != '')
       ]
       
       if len(df_filtered) == 0:
           logging.warning(f"No rows found for {cdr_region} in peptide {peptide}")
           continue
       
       # Store embeddings and indices for vurrent CDR region
       region_embeddings = []
       region_raw_indexes = []
       
       # Extract embeddings for each sequence
       for _, row in df_filtered.iterrows():
           try:
               # Find sequence in embeddings using raw-index
               idx = raw_indexes.tolist().index(row['raw_index'])
               full_embedding = embeddings[idx]
               
               # Locate CDR sequence in full chain
               chain_seq = row[chain_type]
               start_pos, end_pos = find_cdr_positions(chain_seq, row[cdr_region])
               
               if start_pos is not None:
                   # Extract embedding only for CDR 
                   cdr_embedding = full_embedding[start_pos+1:end_pos+1]
                   logging.info(f"CDR sequence length: {len(row[cdr_region])}, Embedding shape: {cdr_embedding.shape}")
                   
                   region_embeddings.append(cdr_embedding)
                   region_raw_indexes.append(row['raw_index'])
               else:
                   logging.warning(f"Could not find CDR sequence in {chain_type} for raw_index {row['raw_index']}")
                   
           except ValueError:
               logging.warning(f"Raw index {row['raw_index']} not found in embeddings")
       
       # Save extracted embeddings
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
       # Extract embeddings for alpha chain 
       alpha_cdrs = ['CDR1a', 'CDR2a', 'CDR3a']
       extract_embeddings_for_chain(peptide, 'TCRa', alpha_cdrs)
       
       # Extract embeddings for beta chain
       beta_cdrs = ['CDR1b', 'CDR2b', 'CDR3b']
       extract_embeddings_for_chain(peptide, 'TCRb', beta_cdrs)
       
       logging.info(f"Processing completed for peptide {peptide}")
   except Exception as e:
       logging.error(f"An error occurred: {e}")
       raise

if __name__ == "__main__":
   main()