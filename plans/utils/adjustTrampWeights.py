#!/usr/bin/env python3

import sys
import shutil
import argparse
import os

import tramp
import opt4d_bwf

def setWeight(tramp_data, bwf, scale):
    print('Total weight = {}'.format(bwf.totalWeight))
    tramp_data.gigaproton_total = bwf.totalWeight*scale
    for i in range(tramp_data.number_spots):
        tramp_data.spots[i].ngp = bwf.w[i]*scale

def supressZeros(tramp_data):
    temp = []
    for i in range(tramp_data.number_spots):
        if tramp_data.spots[i].ngp >= 0.000000001:
            temp.append(tramp_data.spots[i])
    tramp_data.number_spots = len(temp)
    tramp_data.spots = temp[:]


def adjust(intrampfile, inbwffile, fractions):
    print(os.getcwd())

    print('Reading tramp {}'.format(intrampfile))
    tramp_data = tramp.Tramp(intrampfile)
    print('Reading bwf {}'.format(inbwffile))
    bwf = opt4d_bwf.Bwf(inbwffile)
    if bwf.number_spots != tramp_data.number_spots:
        print("ERROR! The # of spots in the tramp and bwf files are different!!");
        exit();

    print('Setting new weights')
    setWeight(tramp_data, bwf, 1.0/fractions)
    supressZeros(tramp_data)
    
    shutil.copy2(intrampfile, intrampfile+'_backup')
    
    outfile = intrampfile
    print('Writting {}'.format(outfile))
    f = open(outfile, 'w')
    f.writelines(tramp_data.getHeader())
    for i in range(tramp_data.number_spots):
        spot = tramp_data.spots[i]
        f.write(str.format("{0:.3f}", spot.energy)+'\t'+str.format("{0:.3f}", spot.x)+'\t'+str.format("{0:.3f}", spot.y)+'\t'+str.format("{0:.9f}", spot.ngp)+'\n')



def adjust_main(args):
    adjust(args.tramp, args.weights, args.fractions)

def main(argv):
    parser = argparse.ArgumentParser(description='Replace weights in tramp with optimized weights from Opt4D')
    parser.add_argument('--tramp', required=True,
                                help='Input tramp file to fill')
    parser.add_argument('--weights', required=True,
                                help='Input file with optimized weights')
    parser.add_argument('--fractions', default=30, type=int,
                                help='Number of fractions in which the plan is divided.')
    parser.set_defaults(func=adjust_main)
    args = parser.parse_args()
    args.func(args)

if __name__=="__main__":
    main(sys.argv[1:])


