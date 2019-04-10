import numpy as np
import pylab as plt
from fast_tsne import fast_tsne
import seaborn as sns
from sklearn import preprocessing
sns.set()

DeepAE_Original = np.load("./Compressed_data/GSE69405_Original.npy")
DeepAE_encoded_10 = np.load("./Compressed_data/GSE69405_encoded_10.npy")
DeepAE_decoded_10 = np.load("./Compressed_data/GSE69405_decoded_10.npy")
DeepAE_encoded_25 = np.load("./Compressed_data/GSE69405_encoded_25.npy")
DeepAE_decoded_25 = np.load("./Compressed_data/GSE69405_decoded_25.npy")
DeepAE_encoded_50 = np.load("./Compressed_data/GSE69405_encoded_50.npy")
DeepAE_decoded_50 = np.load("./Compressed_data/GSE69405_decoded_50.npy")

print(DeepAE_Original.shape)
print(DeepAE_encoded_10.shape)
print(DeepAE_decoded_10.shape)
print(DeepAE_encoded_25.shape)
print(DeepAE_decoded_25.shape)
print(DeepAE_encoded_50.shape)
print(DeepAE_decoded_50.shape)

#print(DeepAE_Original[:2, :10])
#print(DeepAE_encoded_10[:2, :10])
#print(DeepAE_decoded_10[:2, :10])
#print(DeepAE_encoded_25[:2, :10])
#print(DeepAE_decoded_25[:2, :10])
#print(DeepAE_encoded_50[:2, :10])
#print(DeepAE_decoded_50[:2, :10])

# col = np.array(['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99',
#                '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a'])

Z1 = fast_tsne(DeepAE_Original, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z1 = preprocessing.minmax_scale(Z1)

Z2 = fast_tsne(DeepAE_encoded_10, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z2 = preprocessing.minmax_scale(Z2)

Z3 = fast_tsne(DeepAE_decoded_10, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z3 = preprocessing.minmax_scale(Z3)

Z4 = fast_tsne(DeepAE_encoded_25, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z4 = preprocessing.minmax_scale(Z4)

Z5 = fast_tsne(DeepAE_decoded_25, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z5 = preprocessing.minmax_scale(Z5)

Z6 = fast_tsne(DeepAE_encoded_50, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z6 = preprocessing.minmax_scale(Z6)

Z7 = fast_tsne(DeepAE_decoded_50, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
Z7 = preprocessing.minmax_scale(Z7)

fig = plt.figure(figsize=(10, 10))
plt.subplot(331)
plt.scatter(Z1[:, 0], Z1[:, 1], s=1)
plt.title('Original_data(GSE69405)')
plt.subplot(332)
plt.scatter(Z2[:, 0], Z2[:, 1], s=1)
plt.title('DeepAE_Encode_10')
plt.subplot(333)
plt.scatter(Z3[:, 0], Z3[:, 1], s=1)
plt.title('DeepAE_Decode_10')
plt.subplot(335)
plt.scatter(Z4[:, 0], Z4[:, 1], s=1)
plt.title('DeepAE_Encode_25')
plt.subplot(336)
plt.scatter(Z5[:, 0], Z5[:, 1], s=1)
plt.title('DeepAE_Decode_25')

plt.subplot(338)
plt.scatter(Z6[:, 0], Z6[:, 1], s=1)
plt.title('DeepAE_Encode_50')
plt.subplot(339)
plt.scatter(Z7[:, 0], Z7[:, 1], s=1)
plt.title('DeepAE_Decode_50')

plt.show()
fig.savefig('./t-SNE_figures/GSE69405_2.png', bbox_inches='tight')

#np.save('Z1.npy', Z1)
#np.save('Z2.npy', Z2)
#np.savetxt("Z1.csv", Z1, delimiter=',')
#np.savetxt("Z2.csv", Z2, delimiter=',')
