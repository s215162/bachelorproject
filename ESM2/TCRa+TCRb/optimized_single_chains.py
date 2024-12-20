import torch
import esm
import pandas as pd
import os
import logging
import sys
import numpy as np

def setup_logging(peptide, sequence_type, paths):
    """Configure logging to both file and console output for tracking ESM processing"""
    log_file = os.path.join(paths['log_dir'], 
                           f"optimized_{sequence_type}_ESM_{peptide}_matrix_processing.log")
    # configuration for file logging
    logging.basicConfig(
        filename=log_file,
        filemode="a",
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )
    # Add console outpot for monitoring
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

def process_sequences(df_filtered, peptide, sequence_type, paths, 
                     device='cuda' if torch.cuda.is_available() else 'cpu'):
    """Generate ESM-2 embeddings for TCR sequences and save as numpy arrays
    
    Uses ESM-2 (650M parameter model) to create sequence embeddings for each TCR.
    Processes sequences in batches to manage memory usage efficiently.
    """
    # Check if there are sequences to process
    if len(df_filtered) == 0:
        logging.error(f"No sequences found for peptide {peptide}")
        return
    
    # Initialise storage for embeddings and their indices
    sequence_representations = []
    raw_indexes = []
    # to ensure a larger batch size for GPU to improve processing speed
    batch_size = 50 if torch.cuda.is_available() else 25

    # Load and prepare ESM-2 model
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device)
    model.eval()  # Set model to evaluation mode
    batch_converter = alphabet.get_batch_converter()

    all_data = [(row["peptide_x"], row[sequence_type]) for _, row in df_filtered.iterrows()]
    all_raw_indexes = df_filtered["raw_index"].tolist()

    # Process sequences in batches to help memory capacity
    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        
        logging.info(f"Processing batch {start//batch_size + 1}/{(len(df_filtered)-1)//batch_size + 1}")
        
        try:
            batch_data = all_data[start:end]
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            
            # Move data to appropriate device (GPU/CPU)
            batch_tokens = batch_tokens.to(device)
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33], return_contacts=True)

            # extract embeddings from layer 33 (final layer)
            token_representations = results["representations"][33]

            # Store embeddings for each sequence, excluding start/end tokens
            for i, tokens_len in enumerate(batch_lens):
                seq_rep = token_representations[i, 1:tokens_len - 1].cpu().numpy()
                sequence_representations.append(seq_rep)
                raw_indexes.append(all_raw_indexes[start + i])

        except Exception as e:
            logging.error(f"Error processing batch {start}-{end}: {e}")
            continue

        # periodically clear GPU memory to prevent OOM errors
        if device == 'cuda' and (start//batch_size) % 10 == 0:
            torch.cuda.empty_cache()

    if not sequence_representations:
        logging.error("No sequence representations were generated")
        return

    # save embeddinfs and indices to npz file
    npz_file = os.path.join(paths['output_dir'], 
                           f"final_esm_embedding_matrix_{sequence_type}_data_{peptide}.npz")
    try:
        np.savez(npz_file, 
                 embeddings=np.array(sequence_representations, dtype=object), 
                 raw_indexes=raw_indexes)
        logging.info(f"Saved {len(sequence_representations)} embeddings to {npz_file}")
    except Exception as e:
        logging.error(f"Failed to save embeddings: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <peptide> <sequence_type>")
        sys.exit(1)

    peptide = sys.argv[1]
    sequence_type = sys.argv[2]  # Can be TCRb, TCRa, or tcr_full
    
    # Configure processing paths - update to user paths
    paths = {
        'data_dir': '',    # sequence data directory
        'output_dir': '',  # output directory
        'log_dir': ''     # log directory
    }
    
    setup_logging(peptide, sequence_type, paths)

    try:
        # load and filter dataset
        logging.info(f"Processing {sequence_type} sequences for peptide {peptide}")
        df = pd.read_csv(os.path.join(paths['data_dir'], 'full_sequence_data.csv'))
        df_filtered = df[df["peptide_x"] == peptide].copy()  
        
        # generate and save embeddings
        process_sequences(df_filtered, peptide, sequence_type, paths)
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()