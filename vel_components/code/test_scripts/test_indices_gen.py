import numpy as np

N = 100
sample_size = 50
seed = 1337
rng = np.random.default_rng(seed)
sample_indices = np.empty(0)
while sample_indices.size != sample_size:
    random_indices = rng.integers(N, size=(sample_size - sample_indices.size))
    sample_indices = np.unique(np.concatenate((sample_indices, random_indices)))
    seed += 1
    print("Try {}, Seed is now {}, # of Indices is {}".format(seed - 1337, seed, sample_indices.size))
    print("Sample Indices:")
    print(sample_indices)
    print()