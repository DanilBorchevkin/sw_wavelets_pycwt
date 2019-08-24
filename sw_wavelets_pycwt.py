import glob
import os
import pycwt
#from __future__ import division
import numpy
from matplotlib import pyplot
import pycwt as wavelet
from pycwt.helpers import find

def get_graph_from_file(in_filepath, out_folder, out_filename):
    # Get data
    # TODO there are differents formats of file
    # TODO implement differents parsers by parameters of function
    p1 = numpy.genfromtxt(in_filepath)

    # TODO fix this shit
    dat = p1

    # TODO tihnk about kwargs for this values
    title = "Demo title"
    label = "Demo label"
    units = "Demo unit"

    # Values for calculations
    # TODO spike about args
    t0 = 0          # start time
    dt = 1          # step of differentiation

    # Time array in years - TODO not in yers. think about it
    N = dat.size
    t = numpy.arange(0, N) * dt + t0

    # Differentiation
    # TODO use enumarete and check odd and evens
    i = 0
    dat_diff_list = list()
    for val in dat:
        if i == 0:
            i += 1
            continue
        
        cur_diff = (dat[i] - dat[i-1]) / dt
        dat_diff_list.append(cur_diff)

        i += 1
    dat_diff = numpy.asarray(dat_diff_list)
    t_diff = t[:299]

    # Detrenting - TODO not neccessary
    # TODO think about neccessary
    dat_norm = dat_diff
    
    # Define parameters for wavelet
    omega_null = 6
    mother = wavelet.Morlet(omega_null)
    s0 = 2 * dt     # starting scale
    dj = 1 / 12     # amount of suboctaves in octave
    J = 7 / dj      # Seven powers of two with dj sub-octaves
    alpha, _,_ = wavelet.ar1(dat_diff)  # Lag-1 correlation for red noise

    # Perform wavelet transform
    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J, mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother)# * std
    # Inverse wavelet transform # TODO we need normalization
    #iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std


    # Plot
    pyplot.close('all')
    ax = pyplot.axes()
    
    #pyplot.ioff()
    #figprops = dict(figsize=(11, 8), dpi=72)
    #fig = pyplot.figure(**figprops)

    ax.plot(t_diff, iwave, '-', linewidth=1, color=[0.5, 0.5, 0.5])
    ax.plot(t_diff, dat_diff)

    # Save graph to file
    # TODO implement
    #plt.savefig('{}/{}.png'.format(out_folder, out_filename))
    # ----------------------------------------------
    # or show the graph
    pyplot.show()

def get_graph_name_from_filepath(filepath):
    result = ""

    # Get only file name
    result = os.path.basename(filepath)

    # Truncate extension
    result = result[: result.rfind(".")]

    # Truncate all before station name
    #result = result[result.find("wevelet_") + len("wevelet_") :]

    # Truncate passage number
    # Some tricky moment we need second from right _
    #third_ = result.rfind("_")
    #second_ = result.rfind("_", 0, third_) 
    #result = result[:second_] + result[third_:]

    return result

def process_all_files_in_folder(in_folder, file_mask, out_folder):
    for filepath in glob.glob(in_folder + "/" + file_mask):
        graph_name = get_graph_name_from_filepath(filepath)
        print(" >>> Process file '{}' with output name '{}'".format(filepath, graph_name))
        get_graph_from_file(filepath, out_folder, graph_name)

def main():
    print("Script is started")
    process_all_files_in_folder("./input", "*.txt", "./output")

if __name__ == "__main__":
    main()