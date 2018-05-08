'''
Plot strong scaling data.
'''

from matplotlib import pyplot as plt

# dim = 30000
# sparsity = 100
# iters = 100

data = {
    "nprocs" : [
        1,
        2,
        4,
        8,
        12,
        16,
        24,
        32
    ],
    "PETSc" : [
        5.286902,   # 1
        2.772392,   # 2
        0.881947,   # 4
        0.375177,   # 8
        0.253327,   # 12
        0.199869,   # 16
        0.155552,   # 24
        0.141784,   # 32
    ],
    "BlockCSR" : [
        4.42098,    # 1
        2.80139,    # 2
        0.934391,   # 4
        0.431943,   # 8
        0.302391,   # 12
        0.204759,   # 16
        0.160694,   # 24
        0.12086,    # 32
    ],
    "SingleCSR" : [
        3.64,       # 1
        180.523,    # 2
        140.28,     # 4
        81.799,     # 8
        78.9919,    # 12
        59.6613,    # 16
        40.7031,    # 24
        30.6574,    # 32
    ],
    "RC" : [
        5.21182,    # 1
        3.36233,    # 2
        2.84909,    # 4
        2.5633,     # 8
        2.59377,    # 12
        2.62786,    # 16
        2.72344,    # 24
        2.74589,    # 32
    ],
}

# compare parallel efficiency to 1-process PETSc
baseline = data['PETSc'][0]

linestyles = {
    "RC" : '--',
    "SingleCSR" : ':',
    "BlockCSR" : '-',
    "PETSc" : '-.',
}

f, ax = plt.subplots()
f.set_size_inches(5,4)

for k in linestyles:
    plt.plot(data['nprocs'], [baseline/(v*p) for v,p in zip(data[k], data['nprocs'])],
             color='0.0', linestyle=linestyles[k], marker='.', label=k)

plt.xlabel("Number of Processors")
plt.ylabel("Parallel Efficiency [Inverse CPU-s, scaled to 1-process PETSc]")

ymin, ymax = plt.ylim()
plt.ylim(0,ymax)

plt.legend()

plt.show()
