intersect = np.intersect1d(a,b)
asort = argsort(a)
a_ind = np.searchsorted(a,intersect,sorter = asort)
bsort = argsort(b)
b_ind = np.searchsorted(b,intersect,sorter = bsort)
inda = asort[a_ind]
indb = bsort[b_ind]
