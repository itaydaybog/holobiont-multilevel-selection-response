import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import mls_general_code_original as mlsg
import numpy.lib.recfunctions as rf
from pathlib import Path
import os
from scipy import interpolate
from mpl_toolkits.mplot3d import axes3d

data_folder = Path("Step_Data/")
fig_Folder = Path("Step_Figures/")
figureName = 'figure3_fitness_step_function_full.pdf'

sigma_vec = [0.02, 0.1]
tauVRange = (-2, 2)
tauHRange = (-2, 4)

interpolate_2d_graph = True

"""
# SET figure settings
"""
wFig = 8.7
hFig = 4.5
font = {'family': 'Helvetica',
        'weight': 'light',
        'size': 6}

axes = {'linewidth': 0.5,
        'titlesize': 7,
        'labelsize': 6,
        'labelpad': 2,
        'spines.top': False,
        'spines.right': False,
        }

ticks = {'major.width': 0.5,
         'direction': 'in',
         'major.size': 2,
         'labelsize': 6,
         'major.pad': 2}

legend = {'fontsize': 6,
          'handlelength': 1.5,
          'handletextpad': 0.5,
          'labelspacing': 0.2}

figure = {'dpi': 300}
savefigure = {'dpi': 300,
              'transparent': True}

mpl.style.use('seaborn-ticks')
mpl.rc('font', **font)
mpl.rc('axes', **axes)
mpl.rc('xtick', **ticks)
mpl.rc('ytick', **ticks)
mpl.rc('legend', **legend)
mpl.rc('figure', **figure)
mpl.rc('savefig', **savefigure)


# Process data, calc timescale and other variables needed for plotting
def process_data(statData):
    # calculate heritability time
    tauHer = mlsg.calc_tauHer_numeric(
        statData['n0'], statData['mig'])
    tauVar = mlsg.calc_tauV(statData['cost'])
    tauHerRel = tauHer/statData['TAU_H']
    tauVar_rel = tauVar/statData['TAU_H']
    sigma_cat = mlsg.make_categorial(statData['sigmaBirth'])
    BH_cat = mlsg.make_categorial(statData['B_H'])
    dataToStore = (tauHer, tauVar, tauHerRel, tauVar_rel, sigma_cat, BH_cat)
    nameToStore = ('tauHer', 'tauVar', 'tauHer_rel',
                   'tauVar_rel', 'sigma_cat', 'BH_cat')

    statData = rf.append_fields(
        statData, nameToStore, dataToStore, usemask=False)

    return statData


# get subset of data to plot
def select_data(data1D, BHidx, sigmaidx):
    curSigma = data1D['sigma_cat'] == sigmaidx
    curBH = data1D['BH_cat'] == BHidx
    # remove nan and inf
    isFinite = np.logical_and.reduce(
        np.isfinite((data1D['tauVar_rel'], data1D['tauHer_rel'],
                     data1D['F_mav'], curSigma, curBH)))
    currSubset = np.logical_and.reduce((curSigma, curBH, isFinite))
    # extract data and log transform x,y
    x = np.log10(data1D['tauVar_rel'][currSubset])

    transMode = data1D['n0']/data1D['mig']
    y = np.log10(transMode[currSubset])
    z = data1D['F_mav'][currSubset]
    return (x, y, z)


def plot_3D(ax, data1D, sigmaIndex):
    BHidx = 1
    x, y, z = select_data(data1D, BHidx, sigmaIndex)

    ax.plot_trisurf(x, y, z, linewidth=0, antialiased=False, cmap=mpl.cm.get_cmap('plasma'))

    steps = (3, 4, 3)
    fRange = (0, 1)

    ax.set_xlim(tauVRange)
    ax.set_ylim(tauHRange)
    ax.set_zlim(fRange)
    ax.set_xticks(np.linspace(*tauVRange, steps[0]))
    ax.set_yticks(np.linspace(*tauHRange, steps[1]))
    ax.set_zticks(np.linspace(*fRange, steps[2]))

    # set labels
    ax.set_xlabel('$log_{10} \\frac{\\tau_{Var}}{\\tau_H}$')
    ax.set_ylabel('$log_{10} \\frac{n_0/k}{\\theta/\\beta}$')
    ax.set_zlabel('$\\langle f \\rangle$')
    ax.yaxis.labelpad = -10
    ax.xaxis.labelpad = -10
    ax.zaxis.labelpad = -10
    ax.tick_params(axis='z', which='major', pad=0)
    ax.tick_params(axis='both', which='major', pad=-5)
    ax.view_init(30, -115)

    return None


# bin scatter data to show as heatmap
def bin_2Ddata(currXData, currYData, currZData, xbins, ybins):
    """[Bins x,y data into 2d bins]
    Arguments:
            currXData {np vector} -- xData to bin
            currYData {np vector} -- yData to bin
            currZData {np vector} -- zData to bin
            xbins {np vector} -- xBins to use
            ybins {np vector} -- yBins to use
    """
    # init output
    nX = xbins.size
    nY = ybins.size
    binnedData = np.full((nY, nX), np.nan)
    # loop over bins and calc mean
    for xx in range(nX - 1):
        for yy in range(nY - 1):
            # find data in bin
            inXBin = np.logical_and(
                (currXData >= xbins[xx]), (currXData < xbins[xx+1]))
            inYBin = np.logical_and(
                (currYData >= ybins[yy]), (currYData < ybins[yy+1]))
            inBin = np.logical_and(inXBin, inYBin)
            zInBin = currZData[inBin]
            # calc mean over bine
            binnedData[yy, xx] = np.nanmean(zInBin)
    return(binnedData)


# make heatmap
def plot_heatmap(fig, ax, data1D, sigmaIndex):
    if interpolate_2d_graph:
        lcl_tauVRange = (tauVRange[0], tauVRange[1])
        lcl_tauHRange = (tauHRange[0], tauHRange[1])
    else:
        lcl_tauVRange = (tauVRange[0], tauVRange[1] + 1)
        lcl_tauHRange = (tauHRange[0], tauHRange[1] + 1)
    BHidx = 1
    xStep = 0.25
    yStep = 0.5
    xbins = np.linspace(*lcl_tauVRange, int(
        np.ceil((lcl_tauVRange[1]-lcl_tauVRange[0])/xStep))+1)
    ybins = np.linspace(*lcl_tauHRange, int(
        np.ceil((lcl_tauHRange[1] - lcl_tauHRange[0]) / yStep)) + 1)

    # get data with selection
    xS, yS, zS = select_data(data1D, BHidx, sigmaIndex)
    binnedDataS = bin_2Ddata(xS, yS, zS, xbins, ybins)
    if interpolate_2d_graph:
        x = np.arange(0, binnedDataS.shape[1])
        y = np.arange(0, binnedDataS.shape[0])
        # mask empty values
        binnedDataS = np.ma.masked_invalid(binnedDataS)
        xx, yy = np.meshgrid(x, y)
        # get only the valid values
        x1 = xx[~binnedDataS.mask]
        y1 = yy[~binnedDataS.mask]
        binnedDataS = binnedDataS[~binnedDataS.mask]
        binnedDataS = interpolate.griddata((x1, y1), binnedDataS.ravel(), (xx, yy), method='nearest')
        ax.set_xlim(tauVRange[0], tauVRange[1])
        ax.set_ylim(tauHRange[0], tauHRange[1])
    else:
        ax.set_xlim(tauVRange[0], tauVRange[1] + 0.25)
        ax.set_ylim(tauHRange[0], tauHRange[1] + 0.5)
    im = ax.pcolormesh(xbins, ybins, binnedDataS, cmap='plasma', vmin=0, vmax=1)
    #cb = fig.colorbar(im, ax=ax)
    name = '$\\langle f \\rangle$'
    fig.colorbar(im, ax=ax, orientation='vertical',
                 label=name,
                 ticks=[0, 0.5, 1])

    steps = (3, 4)

    ax.set_xticks(np.linspace(*tauVRange, steps[0]))
    ax.set_yticks(np.linspace(*tauHRange, steps[1]))

    # set labels
    ax.set_xlabel('$log_{10} \\frac{\\tau_{Var}}{\\tau_H}$')
    ax.set_ylabel('$log_{10} \\frac{n_0/k}{\\theta/\\beta}$')
    # plt.show()
    return None


def draw_graph():
    files = os.listdir("Step_Data/")
    data_list = []
    for file in files:
        data_file = np.load(data_folder / file, allow_pickle=True)
        data_list.append(data_file['statData'])
    data1D = np.concatenate(data_list)
    data1D = process_data(data1D)

    fig = plt.figure()
    mlsg.set_fig_size_cm(fig, wFig, hFig)

    # setup manual axis for subplots
    bm = 0.15
    tm = 0.06
    cm = 0.13
    hf = 0.65
    h = [hf * (1 - bm - cm - tm), (1 - hf) * (1 - bm - cm - tm)]
    lm = 0.05
    rm = 0.05
    cmh = 0.1
    w = (1 - cmh - rm - lm) / 2
    wf = 0.85
    wo = (1 - wf) * w
    left = [lm, lm + cmh + w]
    bot = [bm + h[1] + cm, bm]

    # plot for different sampling variance
    for ss in range(2):
        # plot scatter
        ax = fig.add_axes([left[ss], bot[0], w, h[0]], projection='3d')
        plot_3D(ax, data1D, ss)
        ax.set_title('$\\sigma={}$'.format(sigma_vec[ss]), fontsize=6)
        # plot heatmap
        ax = fig.add_axes([left[ss] + wo, bot[1], wf * w, h[1]])
        plot_heatmap(fig, ax, data1D, ss)

    # plt.tight_layout(pad=1, h_pad=2.5, w_pad=1.5)
    fig.savefig(fig_Folder / figureName,
                format="pdf", transparent=True)

    return None


if __name__ == '__main__':
    draw_graph()
