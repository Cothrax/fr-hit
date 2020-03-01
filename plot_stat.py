import os
import sys
import time
import struct
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def read_stat(fn):
    size = os.path.getsize(fn)
    print(fn, size)
    
    arr = np.zeros(size, dtype=np.int8)
    with open(fn, 'rb') as f:
        f.readinto(arr)

    counter = Counter(arr)
    ret = np.zeros(max(arr) + 1)
    for k, v in counter.items():
        ret[k] = v
    return ret


def plot(matches_fn, hits_fn, fig_fn=None):
    if not fig_fn:
        fig_fn = 'fig%s.jpg' % int(time.time())
    
    matches = read_stat(matches_fn)
    hits = read_stat(hits_fn)
    matches_std= matches / np.max(matches)
    hits_std = hits / np.max(hits)
    width = 0.3

    plt.bar(np.arange(len(matches_std)) - width, matches_std, width=width, align='edge')
    plt.bar(np.arange(len(hits_std)) + width, hits_std, width=width, align='edge')
    plt.plot(np.cumsum(matches) / np.sum(matches))
    plt.plot(np.cumsum(hits) / np.sum(hits))

    plt.savefig(fig_fn)
    plt.show()


    
if __name__ == '__main__':
    print('in', os.path.pardir)
    plot(*(sys.argv[1:]))
