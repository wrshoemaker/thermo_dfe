from __future__ import division
import os
import pandas as pd
import  matplotlib.pyplot as plt

mydir = os.path.expanduser("~/GitHub/thermo_dfe/")

def make_hist():
    IN_path = mydir + 'data/deltas.txt'
    IN = pd.read_csv(IN_path, sep = '\t')
    s = IN.Selection.values
    fig = plt.figure()
    plt.hist(s, bins=100, alpha = 0.8,  normed = True)
    #plt.title(sample + ' dist. of coverage')
    plt.xlabel('Selection coefficient', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.xlim([min(s), max(s)])
    fig.tight_layout()
    out_plot = mydir + 'figures/hist.png'
    fig.savefig(out_plot, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def make_hist_zoom():
    IN_path = mydir + 'data/deltas.txt'
    IN = pd.read_csv(IN_path, sep = '\t')
    s = IN.Selection.values
    s = [s_i for s_i in s if s_i > 0]
    fig = plt.figure()
    plt.hist(s, bins=100, alpha = 0.8,  normed = True)
    #plt.title(sample + ' dist. of coverage')
    plt.xlabel('Selection coefficient', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.xlim([0, max(s)])
    print max(s)
    fig.tight_layout()
    out_plot = mydir + 'figures/hist_zoom.png'
    fig.savefig(out_plot, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

make_hist()
make_hist_zoom()
