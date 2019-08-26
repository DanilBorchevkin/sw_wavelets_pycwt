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

    title = 'NINO3 Sea Surface Temperature'
    label = 'NINO3 SST'
    units = 'degC'

    # Values for calculations
    # TODO spike about args
    t0 = 12           # start time
    dt = 0.5          # step of differentiation - in minutes

    N = dat.size
    t = numpy.arange(0, N) * dt + t0

    p = numpy.polyfit(t - t0, dat, 1)
    dat_notrend = dat - numpy.polyval(p, t - t0)
    std = dat_notrend.std()  # Standard deviation
    var = std ** 2  # Variance
    dat_norm = dat_notrend / std  # Normalized dataset

    mother = wavelet.Morlet(6)
    s0 = 2 * dt  # Starting scale, in this case 2 * 0.25 years = 6 months
    dj = 1 / 12  # Twelve sub-octaves per octaves
    J = 7 / dj  # Seven powers of two with dj sub-octaves
    alpha, _, _ = wavelet.ar1(dat)  # Lag-1 autocorrelation for red noise

    wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(dat_norm, dt, dj, s0, J, mother)
    iwave = wavelet.icwt(wave, scales, dt, dj, mother) * std

    power = (numpy.abs(wave)) ** 2
    fft_power = numpy.abs(fft) ** 2
    period = 1 / freqs

    power /= scales[:, None]

    signif, fft_theor = wavelet.significance(1.0, dt, scales, 0, alpha, significance_level=0.95, wavelet=mother)
    sig95 = numpy.ones([1, N]) * signif[:, None]
    sig95 = power / sig95

    glbl_power = power.mean(axis=1)
    dof = N - scales  # Correction for padding at edges
    glbl_signif, tmp = wavelet.significance(var, dt, scales, 1, alpha, significance_level=0.95, dof=dof, wavelet=mother)

    sel = find((period >= 2) & (period < 8))
    Cdelta = mother.cdelta
    scale_avg = (scales * numpy.ones((N, 1))).transpose()
    scale_avg = power / scale_avg  # As in Torrence and Compo (1998) equation 24
    scale_avg = var * dj * dt / Cdelta * scale_avg[sel, :].sum(axis=0)
    scale_avg_signif, tmp = wavelet.significance(var, dt, scales, 2, alpha, significance_level=0.95, dof=[scales[sel[0]], scales[sel[-1]]], wavelet=mother)

    # Prepare the figure
    pyplot.close('all')
    #pyplot.ioff()
    figprops = dict(dpi=144)
    fig = pyplot.figure(**figprops)

    # Second sub-plot, the normalized wavelet power spectrum and significance
    # level contour lines and cone of influece hatched area. Note that period
    # scale is logarithmic.
    bx = pyplot.axes([0.1, 0.37, 0.65, 0.28])
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    bx.contourf(t, numpy.log2(period), numpy.log2(power), numpy.log2(levels),
                extend='both', cmap=pyplot.cm.viridis)
    extent = [t.min(), t.max(), 0, max(period)]
    bx.contour(t, numpy.log2(period), sig95, [-99, 1], colors='k', linewidths=2,
            extent=extent)
    bx.set_title('{} Wavelet Power Spectrum ({})'.format(label, mother.name))
    bx.set_ylabel('Period (minutes)')
    #
    #Yticks = 2 ** numpy.arange(numpy.ceil(numpy.log2(period.min())),
    #                        numpy.ceil(numpy.log2(period.max())))
    #bx.set_yticks(numpy.log2(Yticks))
    #bx.set_yticklabels(Yticks)
    #bx.set_yticks(Yticks)
    #bx.set_yticks(period)

    # Save graph to file
    # TODO implement
    #pyplot.savefig('{}/{}.png'.format(out_folder, out_filename))
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