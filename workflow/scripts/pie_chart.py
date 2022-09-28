#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
import numpy as np

def pie_rbp(ax1):
    values = np.array([31,91])
    labels = ["TPM < 10",r"TPM$\geq$10"]
    explode = [0.01,0.01]
    ax1.pie(values,explode=explode, colors = ['lightgrey','#A6BDDB'], textprops={'fontsize': 24},
        autopct=lambda x: '{:.0f}'.format(x*values.sum()/100))
    patches, texts= ax1.pie(values, colors=['lightgrey','#A6BDDB'],explode=explode)
    ax1.legend(patches, labels, loc=(0.22,-0.2),fontsize=18)
    ax1.set_title('ENCODE RBPs',fontsize=24,y=0.97)
    return ax1

def bar_of_pie():
    # make figure and assign axis objects
    fig, (ax3, ax1,ax2) = plt.subplots(1, 3, figsize=(15, 7))
    #fig.suptitle('Filter Nodes by TPM',fontweight="bold")
    fig.subplots_adjust(wspace=0)

    # pie chart parameters
    y = np.array([381,185])
    labels = ["TPM < 10",r"TPM$\geq$10"]
    explode = [0.01,0.01]
    # rotate so that first wedge is split by the x-axis
    wedges, *_ = ax1.pie(y, startangle=60,explode=explode,colors = ['lightgrey','#FC9272'],autopct=lambda x: '{:.0f}'.format(x*y.sum()/100),textprops={'fontsize': 24})
    ax1.set_title('Box C/D snoRNAs',fontsize=24,y=0.97)
    ax1.legend(wedges, labels, loc=(0.22,-0.2),fontsize=18)

    # bar chart parameters
    age_ratios = [155,30]
    age_labels = ['#copies < 10','#copies$\geq$10']
    #colors = ['lightgrey','salmon']
    bottom = 1
    width = .15

    # Adding from the top matches the legend.
    for j, (height, label) in enumerate(reversed([*zip(age_ratios, age_labels)])):
        bottom -= height
        bc = ax2.bar(0, height, width, bottom=bottom, color=['salmon'], label=label,
                    alpha=0.1 + 0.25 * j)
        ax2.bar_label(bc, labels=[f"{height}"], label_type='center',fontsize=22)

    ax2.set_title('# of snoRNA copies',fontsize=18,y=0.95)
    ax2.legend(loc=(0.63,0),fontsize=13)
    ax2.axis('off')
    ax2.set_xlim(- 2.5 * width, 2.5 * width)

    # use ConnectionPatch to draw lines between the two plots
    theta1, theta2 = wedges[0].theta1, wedges[0].theta2
    center, r = wedges[0].center, wedges[0].r
    bar_height = -185#sum(age_ratios)
    #"""
    # draw top connecting line
    x = r * np.cos(np.pi / 180 * theta2) + center[0]
    y = r * np.sin(np.pi / 180 * theta2) + center[1]
    con = ConnectionPatch(xyA=(-width / 2, bar_height), coordsA=ax2.transData,
                        xyB=(x, y), coordsB=ax1.transData)
    con.set_color([0, 0, 0])
    con.set_linewidth(0.5)
    ax2.add_artist(con)

    # draw bottom connecting line
    x = r * np.cos(np.pi / 180 * theta1) + center[0]
    y = r * np.sin(np.pi / 180 * theta1) + center[1]
    con = ConnectionPatch(xyA=(-width / 2, 0), coordsA=ax2.transData,
                        xyB=(x, y), coordsB=ax1.transData)
    con.set_color([0, 0, 0])
    ax2.add_artist(con)
    con.set_linewidth(0.5)
    #"""
    pie_rbp(ax3)
    #plt.figtext(0.47,0.75,"Total = 566",fontsize=14)
    #plt.figtext(0.21,0.75,"Total = 122",fontsize=14)
    #plt.figtext(0.5,0.05,"* cell lines: HCT116, MCF7, PC3, SKOV_frg, TOV112D",fontsize=10)
    plt.savefig('../../results/node_bar_pie.png')
    return

def main():
    """
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(9,4.5))
    fig.suptitle('Filter Nodes by TPM',fontweight="bold")
    # snnoRNA
    #y = np.array([185,381])
    y = np.array([381,185])
    #labels = [r"TPM$\geq$10 in at least"+"\n"+"one cell line","TPM < 10"+"\n"+"in all cell lines"]
    #labels = ["TPM < 10"+"\n"+"in all cell lines",r"TPM$\geq$10 in at least"+"\n"+"one cell line"]
    labels = ["TPM < 10 in all cell lines",r"TPM$\geq$10 in at least one cell line"]
    explode = [0.01,0.01]

    #ax1.pie(y, labels = labels, explode = explode,colors = ['lightgrey','#FC9272'], startangle=-120,
    #    autopct=lambda x: '{:.0f}'.format(x*y.sum()/100))
    ax1.pie(y, explode = explode,colors = ['lightgrey','#FC9272'], startangle=-120,
        autopct=lambda x: '{:.0f}'.format(x*y.sum()/100))
    ax1.set_title('Box C/D snoRNAs',fontsize=11,y=0.97)

    patches, texts= ax1.pie(y, colors=['lightgrey','#FC9272'], startangle=-120,explode=explode)
    ax1.legend(patches, labels, loc=(0.15,-0.05),fontsize=8)
    pie_rbp(ax2)

    plt.figtext(0.5,0.05,"* cell lines: HCT116, MCF7, PC3, SKOV_frg, TOV112D",fontsize=8)
    plt.figtext(0.258,0.8,"Total = 566",fontsize=9)
    plt.figtext(0.685,0.8,"Total = 122",fontsize=9)
    plt.savefig('../../results/node_by_tpm.png')
    """
    bar_of_pie()


    # rbp = auamarine and silver
if __name__ == '__main__':
    main()