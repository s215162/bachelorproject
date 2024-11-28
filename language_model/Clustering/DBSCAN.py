import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
import seaborn as sns

# Load your dataset
data = pd.read_csv('/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/tcra_tcrb_similarity.csv')

# Normalize the sum_max_similarity
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data[['sum_max_similarity']])

# Run DBSCAN
dbscan = DBSCAN(eps=0.0001, min_samples=5)  # Adjust eps and min_samples as needed
data['cluster'] = dbscan.fit_predict(data_scaled)

# Visualize clusters with binders and non-binders
plt.figure(figsize=(10, 6))
sns.scatterplot(
    x=data['sum_max_similarity'], 
    y=[0] * len(data), 
    hue=data['binder'], 
    style=data['cluster'], 
    palette='Set1'
)
plt.xlabel('Sum Max Similarity')
plt.title(f'DBSCAN Clustering (Binders vs. Non-binders)')

# Save the plot
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/dbscan_clustering.png'
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save
plt.close()  # Close figure

# Evaluate cluster results
num_clusters = len(set(data['cluster'])) - (1 if -1 in data['cluster'] else 0)  # Exclude noise points
print(f"Number of Clusters: {num_clusters}")
noise_points = len(data[data['cluster'] == -1])
print(f"Number of Noise Points: {noise_points}")

# Silhouette Score (only if there are more than 1 clusters and no noise points)
if num_clusters > 1 and -1 not in data['cluster']:
    silhouette_avg = silhouette_score(data_scaled, data['cluster'])
    print(f"Silhouette Score: {silhouette_avg}")
else:
    print("Silhouette Score cannot be calculated for less than 2 clusters or with noise.")

# Save clustered data
data.to_csv('/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/dbscan_clustered_data.csv', index=False)
