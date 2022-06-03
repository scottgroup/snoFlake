#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

def main():

    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,10))

    # snnoRNA
    y = np.array([164,163,1214])
    labels = ["C/D","H/ACA","TPM < 10"]
    explode = [0,0,0.2]

    ax1.pie(y, labels = labels, explode = explode,colors = ['salmon','silver','silver'])
    
    y2 = np.array([21,91])
    labels = ["others",r"TPM$\geq$10"]
    explode = [0,0.2]
    ax2.pie(y2,labels = labels, explode=explode, colors = ['silver','xkcd:turquoise'])
    plt.show() 


    # rbp = auamarine and silver
if __name__ == '__main__':
    main()