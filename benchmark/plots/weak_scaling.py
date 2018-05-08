'''
Plot weak scaling data.
'''

from matplotlib import pyplot as plt

# sparsity 100
# iterations 50
# dim 30000 for p = 1

data = {
    "nprocs" : [
        1,
        4,
        9,
        12,
        16,
        24,
        32,
    ],
    "BlockCSR" : [
        1.97196,    # 1, d = 30K
        2.76628,    # 4, d = 60K
        3.64278,    # 9, d = 90K
        3.68643,    # 12, d = 103,923
        4.01226,    # 16, d = 120K
        4.02458,    # 24, d = 146,969
        4.02916,    # 32, d = 169,706
    ],
    "RC" : [
        2.62366,    # 1, d = 30K
        4.16499,    # 4, d = 60K
        6.38339,    # 9, d = 90K
        6.64633,    # 12, d = 103,923
        7.46794,    # 16, d = 120K
        9.29513,    # 24, d = 146,969
        10.2702,    # 32, d = 169,706
    ],
    "PETSc" : [
        2.754051,   # 1, d = 30K
        2.756393,   # 4, d = 60K
        3.425108,   # 9, d = 90K
        3.388548,    # 12, d = 103,923
        3.588989,    # 16, d = 120K
        3.616102,    # 24, d = 146,969
        3.535188,    # 32, d = 169,706
    ],
}

linestyles = {
    "RC" : '--',
    "BlockCSR" : '-',
    "PETSc" : '-.',
}

f, ax = plt.subplots()
f.set_size_inches(5,4)

for k in linestyles:
    plt.plot(data['nprocs'], data[k], color='0.0', linestyle=linestyles[k], marker='.', label=k)

plt.xlabel("Number of Processors")
plt.ylabel("Runtime [s]")

ymin, ymax = plt.ylim()
plt.ylim(0,ymax)

plt.legend()

plt.show()
