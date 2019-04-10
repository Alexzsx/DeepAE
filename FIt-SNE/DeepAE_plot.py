import numpy as np
import pylab as plt
from fast_tsne import fast_tsne
import seaborn as sns
from sklearn import preprocessing
sns.set()

DeepAE_Original = np.load("./Compressed_data/DeepAEGSE60361_Original.npy")
DeepAE_encoded_10 = np.load("./Compressed_data/DeepAEGSE60361_encoded_10.npy")
DeepAE_decoded_10 = np.load("./Compressed_data/DeepAEGSE60361_decoded_10.npy")
DeepAE_encoded_25 = np.load("./Compressed_data/DeepAEGSE60361_encoded_25.npy")
DeepAE_decoded_25 = np.load("./Compressed_data/DeepAEGSE60361_decoded_25.npy")
DeepAE_encoded_50 = np.load("./Compressed_data/DeepAEGSE60361_encoded_50.npy")
DeepAE_decoded_50 = np.load("./Compressed_data/DeepAEGSE60361_decoded_50.npy")

print(DeepAE_Original.shape)
print(DeepAE_encoded_10.shape)
print(DeepAE_decoded_10.shape)
print(DeepAE_encoded_25.shape)
print(DeepAE_decoded_25.shape)
print(DeepAE_encoded_50.shape)
print(DeepAE_decoded_50.shape)


Z1 = fast_tsne(DeepAE_Original, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
#Z1 = Z1[Z1[:, 1] > -5, :]
#Z1 = Z1[Z1[:, 1] < 5, :]
#Z1 = Z1[Z1[:, 0] > -5, :]
#Z1 = Z1[Z1[:, 0] < 5, :]
Z1 = preprocessing.minmax_scale(Z1)

Z2 = fast_tsne(DeepAE_encoded_10, perplexity=50, stop_early_exag_iter=1000, seed=42, nbody_algo='Barnes-Hut')
#Z2 = Z2[Z2[:, 1] > -5, :]
#Z2 = Z2[Z2[:, 1] < 5, :]
#Z2 = Z2[Z2[:, 0] > -5, :]
#Z2 = Z2[Z2[:, 0] < 5, :]
Z2 = preprocessing.minmax_scale(Z2)

fig = plt.figure(figsize=(6, 3))
plt.subplot(121)
plt.scatter(Z1[:, 0], Z1[:, 1], s=1)
plt.title('Original_data(GSE60361)')
#plt.xlim(-1, 1)
#plt.ylim(-0.7, 0.9)

plt.subplot(122)
plt.scatter(Z2[:, 0], Z2[:, 1], s=1)
plt.title('DeepAE(10 measurements)')
#plt.xlim(-1, 1)
#plt.ylim(-0.7, 0.9)

plt.show()
fig.savefig('./t-SNE_figures/GSE60361.pdf', dpi=100, bbox_inches='tight')


