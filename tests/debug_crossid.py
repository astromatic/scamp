#!/usr/bin/env python

from astropy.io import fits
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fname = "../scamp_backup/full_1.cat"
hdu = fits.open(fname)
data = hdu[2].data
data_col = hdu[2].columns
print data_col.names
print len(data)


print "hello"
def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
#for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
#    xs = randrange(n, 23, 32)
#    ys = randrange(n, 0, 100)
#    zs = randrange(n, zlow, zhigh)
#    ax.scatter(xs, ys, zs, c=c, marker=m)

d = dict()
for i, val in enumerate(data):
    snum = val['SOURCE_NUMBER']
    if snum not in d:
        d[snum] = list()
    d[snum].append(val)

g = 0
for key, val in d.iteritems():
    print g
    ss = len(val)
    epoch = 0
    c = 'r'
    m = '^'
    if len(val) < 5:
        continue
    g += 1
    for j, val2 in enumerate(val):
        try:
            val3 = val[j+1]
        except:
            ax.scatter(val2['X_IMAGE'], epoch, val2['Y_IMAGE'], c='r')
            break
        ax.plot([val2['X_IMAGE'], val3['X_IMAGE']], [epoch, epoch+1], [val2['Y_IMAGE'], val3['Y_IMAGE']], color='g')
        ax.scatter(val2['X_IMAGE'], epoch, val2['Y_IMAGE'], c='r')
        epoch += 1

    if g > 100:
        break

ax.set_xlabel('X_IMAGE')
ax.set_ylabel('epoch')
ax.set_zlabel('Y_IMAGE')

plt.show()

