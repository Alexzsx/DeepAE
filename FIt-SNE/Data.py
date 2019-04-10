import tensorflow as tf
import numpy as np
import random
from scipy.stats import pearsonr
from scipy.spatial import distance
import seaborn as sns
import matplotlib.pyplot as plt


def Euclidean_dist(A, B):
    C = A - B
    return sum(map(sum, C * C)) ** 0.5


def MAE(A, B):  ## Mean Absolute Error
    C = A - B
    return sum(map(sum, C * C)) / (C.shape[0] * C.shape[1])


def random_split_train_test(X0, training_dictionary_fraction, seed, dictionary_size=0.5, biased_training=0.):
    training_dictionary_size = max(int(training_dictionary_fraction * X0.shape[1]), 5)
    if dictionary_size < 1:
        dictionary_size = dictionary_size * training_dictionary_size
    dictionary_size = int(dictionary_size)
    xi = np.zeros(X0.shape[1], dtype=np.bool)
    if biased_training > 0:
        np.random.seed(seed)
        i = np.random.randint(len(xi))
        dist = distance.cdist([X0[:, i]], X0.T, 'correlation')[0]
        didx = np.argsort(dist)[1:int(biased_training * training_dictionary_size) + 1]
    else:
        didx = []
    xi[didx] = True
    if biased_training < 1:
        remaining_idx = np.setdiff1d(range(len(xi)), didx)
        np.random.seed(seed)
        xi[np.random.choice(remaining_idx, training_dictionary_size - xi.sum(), replace=False)] = True
    xa = X0[:, xi]
    xb = X0[:, np.invert(xi)]
    return xa, xb, xi


def compare_results(A, B):
    results = list((1 - distance.correlation(A.flatten(), B.flatten())))
    results += list(Euclidean_dist(A, B))
    results += list(MAE(A, B))
    return results


seed_all = {"GSE71858": [272, 781, 692, 219, 292],  #
            "GSE60361": [283, 446, 562, 114, 739],  #
            "GSE62270": [629, 685, 953, 595, 378],  #
            "GSE48968": [623, 19, 621, 802, 557],  #
            "GSE52529": [550, 939, 76, 260, 328],  #
            "GSE77564": [475, 649, 316, 639, 741],
            "GSE78779": [152, 866, 808, 796, 184],  #
            "GSE10247": [702, 217, 944, 338, 701],  #
            "GSE69405": [317, 470, 798, 283, 695],
            "GSE45235": [282, 713, 521, 717, 517],  #
            "GSE25038": [480, 402, 413, 64, 574],
            "mass_cytomatry": [943, 800, 175, 486, 749]}

tf.set_random_seed(1)

# Hyper Parameters
LR = 0.0001  # learning rate
Dropout_rate = 0.5
# GSE Data
data_path = "./Original_data/GSE69405.npy"
X = np.load(data_path)
training_dictionary_fraction = 0.05
genes, samples = X.shape
seeds = seed_all['GSE69405']
############################# Define architectures ##################################
# tf placeholder
tf_x = tf.placeholder(tf.float32, [None, genes])  # value in the range of (0, 1)

# encoder
# Dn0 = tf.layers.dropout(tf_x, rate=Dropout_rate, training=True)
en0 = tf.layers.dense(tf_x, 1280, tf.nn.leaky_relu)
en1 = tf.layers.dense(en0, 640, tf.nn.leaky_relu)
en2 = tf.layers.dense(en1, 256, tf.nn.leaky_relu)

encoded = tf.layers.dense(en2, 50)

# decoder
de0 = tf.layers.dense(encoded, 256, tf.nn.leaky_relu)
de1 = tf.layers.dense(de0, 640, tf.nn.leaky_relu)
de2 = tf.layers.dense(de1, 1280, tf.nn.leaky_relu)

decoded = tf.layers.dense(de2, genes, tf.nn.leaky_relu)

loss = tf.losses.mean_squared_error(labels=tf_x, predictions=decoded)
train = tf.train.AdamOptimizer(LR).minimize(loss)

############################# Running ##################################
Results = {}
# seeds = random.sample(range(0, 1000), 5)
# seeds = [283, 446, 562, 114, 739]

print(seeds)
for i in range(1):

    sess = tf.Session()
    sess.run(tf.global_variables_initializer())

    X_train, X_test, xi = random_split_train_test(X, training_dictionary_fraction, seed=seeds[i])
    #print(xi)
    #np.savetxt("GSE60361_Xi.csv", xi, delimiter=',')
    print(X.shape)  #
    print(X_train.shape)  #
    print(X_test.shape)  #
    print(X_train[0, 0:10])
    X_train = np.transpose(X_train)
    X_test = np.transpose(X_test)

    for step in range(500):
        b_x = X_train
        _, encoded_, decoded_, loss_ = sess.run([train, encoded, decoded, loss], {tf_x: b_x})

        if step % 100 == 0:
            # print('------------------Step: %d' % step + '---------------')
            # print('train loss: %.4f' % loss_)
            # plotting decoded image (second row)
            decoded_data_train = sess.run(decoded, {tf_x: b_x})
            # train_p = (1 - distance.correlation(X_train.flatten(), decoded_data_train.flatten()))
            train_pp = pearsonr(X_train.flatten(), decoded_data_train.flatten())[0]
            train_ED = Euclidean_dist(X_train, decoded_data_train)
            train_MAE = MAE(X_train, decoded_data_train)
            # print('train Pearson: %.4f' % train_p)
            # print('train Pearson_: %.4f' % train_pp)
            # print('train Euclidean_dist: %e' % train_ED)
            # print('train MAE: %.4f' % train_MAE)

            encod = sess.run(encoded, {tf_x: b_x})
            # print(encod.shape)
            # print('------------------Test---------------')
            decoded_data_testing = sess.run(decoded, {tf_x: X_test})
            encoded_data = sess.run(encoded, {tf_x: X_test})
            # test_p = (1 - distance.correlation(X_test.flatten(), decoded_data.flatten()))
            test_pp = pearsonr(X_test.flatten(), decoded_data_testing.flatten())[0]
            test_ED = Euclidean_dist(X_test, decoded_data_testing)
            test_MAE = MAE(X_test, decoded_data_testing)
            # print('test Pearson: %.4f' % test_p)
            # print('test Pearson_: %.4f' % test_pp)
            # print('test Euclidean_dist: %e' % test_ED)
            # print('test MAE: %.4f' % test_MAE)
            # print('----------------------------------------')
    #       Result = compare_results(X_test, decoded_data)
    #       print(Result)
    decoded_data_testing = sess.run(decoded, {tf_x: X_test})

    print(decoded_data_testing.shape)

    result_train = 'DeepAE4 (training)_' + str(i)
    result_test = 'DeepAE4 (testing )_' + str(i)
    Results[result_train] = [train_pp, train_ED, train_MAE]
    Results[result_test] = [test_pp, test_ED, test_MAE]
    print('----------------End Iteration: %d' % i + '------------------------')

print(data_path)
for k, v in sorted(Results.items()):
    print('\t'.join([k] + [str(x) for x in v]))

# np.save('DeepAEGSE60361_train.npy', X_train)
#np.save('./Compressed_data/GSE69405_Original.npy', X_test)
np.save('./Compressed_data/GSE69405_encoded_50.npy', encoded_data)
np.save('./Compressed_data/GSE69405_decoded_50.npy', decoded_data_testing)
