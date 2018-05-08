'''
Plot scaling with matrix sparsity.
'''

from matplotlib import pyplot as plt

# sparsity = iterations (possibly scaled)
# dim 10000
# p = 32 nodes, 1 processor per node

data = {
    "sparsity" : [
        1,
        10,
        25,
        100,
        250,
        500,
        750,
        1000,
    ],
    "RC" : [
        0.0151379,  # 1
        0.119727,   # 10
        0.298592,   # 25
        1.12657,    # 100
        2.81762,    # 250
        3.56664,    # 500
        3.5817,     # 750
        3.52355,    # 1000
    ],
    "SingleCSR" : [
        2.67412,    # 1
        3.42723,    # 10
        3.45157,    # 25
        3.45063,    # 100
        3.50171,    # 250
        3.54059,    # 500
        3.607,      # 750
        3.64559,    # 1000
    ],
    "BlockCSR" : [
        0.00930965, # 1
        0.00834375, # 10
        0.0108969,  # 25
        0.0200195,  # 100
        0.0354688,  # 250
        0.065416,   # 500
        0.090024,   # 750
        0.11696,    # 1000
    ],
    "PETSc" : [
        0.00907474, # 1
        0.00826985, # 10
        0.01178220, # 25
        0.02672828, # 100
        0.0583496,  # 250
        0.094636,   # 500
        0.117575,   # 750
        0.168006,   # 1000
    ]
}

linestyles = {
    "RC" : '--',
    "SingleCSR" : ':',
    "BlockCSR" : '-',
    "PETSc" : '-.',
}

f, ax = plt.subplots()
f.set_size_inches(5,4)

for k in linestyles:
    plt.plot([1/s for s in data['sparsity']], data[k],
             color='0.0', linestyle=linestyles[k], marker='.', label=k)

plt.xlabel("Density of nonzeros")
plt.ylabel("Runtime/density [s]")

plt.xscale('log')
plt.yscale('log')

plt.legend()

plt.show()
