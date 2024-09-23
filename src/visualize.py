import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import argparse


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    
    return s == 'True'

parser = argparse.ArgumentParser()
parser.add_argument('--glob', type=boolean_string, default=False)
parser.add_argument('--sa', type=boolean_string, default=False)
parser.add_argument('-show', '--show_layout', type=boolean_string, default=False)
parser.add_argument('-f', '--file', type=str)
args = parser.parse_args()


def str2float(str_list):

    return [float(string) for string in str_list]

def str2int(str_list):

    return [int(string) for string in str_list]


def visualize(file, idx=1):
    glob, sa, show_layout = args.glob, args.sa, args.show_layout
    values = []
    with open(file) as f:
        for line in f:
            values.append(line.split('\n')[0].split(','))
    num_hard_modules = str2int(values[0])[0]
    num_soft_modules = str2int(values[1])[0]
    num_total_modules = num_hard_modules + num_soft_modules
    X = str2float(values[2])
    Y = str2float(values[3])
    W = str2float(values[4])
    H = str2float(values[5])
    Z = str2float(values[6])
    utilization = str2float(values[7])[0]
    bound = str2float(values[8])[0]
    print(values)

    
    label = np.arange(num_total_modules) + 1
    plt.ion()
    fig, ax = plt.subplots()
    for i, txt in enumerate(label):
        if i < num_hard_modules:
            if Z[i] >= 0.9: # Sometimes get 1.01/0.99
                ax.add_patch(Rectangle((X[i], Y[i]), H[i], W[i], color='red'))
                ax.add_patch(Rectangle((X[i], Y[i]), H[i], W[i], color='black', fill=False))
                ax.annotate(text=txt, xy=(X[i], Y[i]), xytext=(X[i]+H[i]/2, Y[i]+W[i]/2))
            else:
                ax.add_patch(Rectangle((X[i], Y[i]), W[i], H[i], color='green'))
                ax.add_patch(Rectangle((X[i], Y[i]), W[i], H[i], color='black', fill=False))
                ax.annotate(text=txt, xy=(X[i], Y[i]), xytext=(X[i]+W[i]/2, Y[i]+H[i]/2))
        else:
            ax.add_patch(Rectangle((X[i], Y[i]), W[i], H[i], color='yellow'))
            ax.add_patch(Rectangle((X[i], Y[i]), W[i], H[i], color='black', fill=False))
            ax.annotate(text=txt, xy=(X[i], Y[i]), xytext=(X[i]+W[i]/2, Y[i]+H[i]/2))
        if sa==True:
            if glob==False:
                    plt.title('Local floorplan for %d-th sub-block: Chip Height = %.4f, Chip Area = %d\nUtilization = %.2f percent' % (idx, bound, bound**2, utilization * 100))
            else:
                    plt.title('Global floorplan for including all sub-blocks: Chip Height = %.4f, Chip Area = %d\nUtilization = %.2f percent' % (bound, bound**2, utilization * 100))
        else:
            plt.title('Direct floorplan: Chip Height = %.4f, Chip Area = %d\nUtilization = %.2f percent' % (bound, bound**2, utilization * 100))

    ax.set_xlim(0, bound)
    ax.set_ylim(0, bound)
    if show_layout:
        plt.show(block=True)
    else:
        plt.close()