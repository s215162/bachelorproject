import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score, confusion_matrix
import seaborn as sns

# Load your dataset
data = pd.read_csv('/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/tcra_tcrb_similarity.csv')


# Normalize the sum_max_similarity for both binders and non-binders
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data[['sum_max_similarity']])

# Run K-means clustering for different values of k
silhouette_scores = []
for k in range(2, 11):
    kmeans = KMeans(n_clusters=k, init='k-means++', random_state=42)
    kmeans.fit(data_scaled)
    silhouette_avg = silhouette_score(data_scaled, kmeans.labels_)
    silhouette_scores.append(silhouette_avg)

# Plot silhouette scores to find the best k
plt.plot(range(2, 11), silhouette_scores)
plt.xlabel('Number of clusters (k)')
plt.ylabel('Silhouette Score')
plt.title('Silhouette Score vs. Number of Clusters')
# Save the plot as an image file
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/silhouetteScore_against_no_of_Clusters.png'  # Change as needed
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save
plt.close()  # Close the figure 

# Perform clustering with the best k
best_k = silhouette_scores.index(max(silhouette_scores)) + 2  # Best k
kmeans = KMeans(n_clusters=best_k, init='k-means++', random_state=42)
data['cluster'] = kmeans.fit_predict(data_scaled)

# Visualize clusters with binders and non-binders
plt.figure(figsize=(10, 6))
sns.scatterplot(x=data['sum_max_similarity'], y=[0] * len(data), hue=data['binder'], style=data['cluster'], palette='Set1')
plt.xlabel('Sum Max Similarity')
plt.title(f'Clustering with {best_k} Clusters (Binders vs. Non-binders)')
# Save the plot as an image file
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/clustering_best_k.png'  # Change as needed
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save
plt.close()  # Close the figure 

# Check how well the clusters segregate binders and non-binders
cluster_purity = []
for cluster in range(best_k):
    cluster_data = data[data['cluster'] == cluster]
    purity = cluster_data['binder'].mean()  # Percentage of binders in the cluster
    cluster_purity.append(purity)

# Display cluster purity (how well binders are separated in each cluster)
print(f"Cluster Purity (Binder Percentage): {cluster_purity}")

# Confusion Matrix for clustering results vs. actual binder labels
conf_matrix = confusion_matrix(data['binder'], data['cluster'])
sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=[f'Cluster {i}' for i in range(best_k)], yticklabels=['Non-binder', 'Binder'])
plt.xlabel('Cluster')
plt.ylabel('Actual Binder Status')
plt.title('Confusion Matrix of Clustering vs. Actual Binder Labels')
# Save the plot as an image file
output_path = '/net/mimer/mnt/tank/projects2/emison/language_model/Clustering/full_similarity/confusion_matrix.png'  # Change as needed
plt.savefig(output_path, dpi=300, bbox_inches='tight')  # High-quality save
plt.close()  # Close figure