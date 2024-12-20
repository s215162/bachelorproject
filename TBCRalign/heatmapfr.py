import scipy
import scipy.cluster.hierarchy as sch
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def create_heatmap(input_path):
   # Load and sort data by peptides
   df = pd.read_csv(input_path)
   print("Data retrieved!")
   
   sorted_df = df.sort_values(by="peptide")
   peptides = sorted_df["peptide"]

   # Find peptide group positions
   peptide_boundaries = [peptides[peptides == pep].index[0] for pep in peptides.unique()]
   peptide_positions = [sorted_df.index.get_loc(boundary) for boundary in peptide_boundaries]

   # Generating labels for peptide the blocks 
   peptide_labels = []
   for i, pep in enumerate(peptides.unique()):
       if i < len(peptide_positions) - 1:
           middle_index = (peptide_positions[i] + peptide_positions[i + 1]) // 2
       else:
           middle_index = (peptide_positions[i] + len(peptides)) // 2
       peptide_labels.append((middle_index, pep))

   # data preperation
   sorted_indices = sorted_df.index.tolist()
   heatmap_data = sorted_df.drop(columns=["peptide", "binder", "partition"])
   heatmap_data.index = sorted_indices
   heatmap_data.columns = sorted_indices
   
   print(heatmap_data.head())

   # create heatmap
   sns.heatmap(heatmap_data)
   plt.xticks(
       [pos for pos, _ in peptide_labels],
       [label for _, label in peptide_labels],
       rotation=90
   )
   plt.yticks(
       [pos for pos, _ in peptide_labels],
       [label for _, label in peptide_labels],
       rotation=0
   )
   plt.title("TBCRalign All CDRs Weighed")
   plt.show()

def main():
   # Add path to your distance matrix file here
   input_path = ""
   create_heatmap(input_path)

if __name__ == "__main__":
   main()