import torch
import esm
import pandas as pd
import os
import logging
import sys
import numpy as np

def setup_logging(peptide, processed_sequence):
    log_file = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/logs/optimized_final_{processed_sequence}_ESM_{peptide}_matrix_cdr_processing.log"
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

def process_sequences(df_filtered, peptide, processed_sequence, device='cuda' if torch.cuda.is_available() else 'cpu'):
    if len(df_filtered) == 0:
        logging.error(f"No sequences found for peptide {peptide}")
        return
    
    sequence_representations = []
    raw_indexes = []
    batch_size = 50 if torch.cuda.is_available() else 25  # Increased batch size for GPU

    # Load ESM-2 model and move to GPU if available
    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    # Pre-process all data at once
    all_data = [(row["peptide_x"], row[processed_sequence]) for _, row in df_filtered.iterrows()]
    all_raw_indexes = df_filtered["raw_index"].tolist()

    # Process in batches
    for start in range(0, len(df_filtered), batch_size):
        end = min(start + batch_size, len(df_filtered))
        
        logging.info(f"Processing batch {start//batch_size + 1}/{(len(df_filtered)-1)//batch_size + 1}")
        
        try:
            # Get current batch
            batch_data = all_data[start:end]
            batch_labels, batch_strs, batch_tokens = batch_converter(batch_data)
            
            # Move batch to GPU/CPU
            batch_tokens = batch_tokens.to(device)
            batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

            with torch.no_grad():
                results = model(batch_tokens, repr_layers=[33], return_contacts=True)

            token_representations = results["representations"][33]

            # Process all sequences in the batch at once
            for i, tokens_len in enumerate(batch_lens):
                seq_rep = token_representations[i, 1:tokens_len - 1].cpu().numpy()
                sequence_representations.append(seq_rep)
                raw_indexes.append(all_raw_indexes[start + i])

        except Exception as e:
            logging.error(f"Error processing batch {start}-{end}: {e}")
            continue

        # Optional: Clear CUDA cache periodically
        if device == 'cuda' and (start//batch_size) % 10 == 0:
            torch.cuda.empty_cache()

    if not sequence_representations:
        logging.error("No sequence representations were generated")
        return

    # Save results
    npz_file_path = f"/net/mimer/mnt/tank/projects2/emison/language_model/final_work/final_esm_embedding_matrix_{processed_sequence}_data_{peptide}.npz"
    try:
        np.savez(npz_file_path, 
                 embeddings=np.array(sequence_representations, dtype=object), 
                 raw_indexes=raw_indexes)
        logging.info(f"Saved {len(sequence_representations)} embeddings to {npz_file_path}")
    except Exception as e:
        logging.error(f"Failed to save embeddings: {e}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <peptide> <sequence_type>")
        sys.exit(1)

    peptide = sys.argv[1]
    # Change depending on needed sequence
    # TCRb  TCRa    tcr_full
    processed_sequence = sys.argv[2]
    setup_logging(peptide, processed_sequence)

    # Define paths
    dataset_path = "/net/mimer/mnt/tank/projects2/emison/language_model/full_sequence_data.csv"
    output_dir = "/net/mimer/mnt/tank/projects2/emison/language_model/final_work/single_chains"
    os.makedirs(output_dir, exist_ok=True)

    try:
        # Load and filter dataset
        logging.info(f"Loading and filtering dataset for peptide {peptide} with sequence {processed_sequence}")
        df = pd.read_csv(dataset_path)
        df_filtered = df[df["peptide_x"] == peptide].copy()  # Added .copy() to avoid SettingWithCopyWarning
        
        # Process sequences
        logging.info(f"Processing dataset for peptide {peptide}")
        process_sequences(df_filtered, peptide, processed_sequence)
        
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

if __name__ == "__main__":
    main()