# --- standard library ---
import os
import math
from datetime import datetime
from itertools import chain

# --- GUI (Tkinter) ---
from tkinter import Tk
from tkinter import *
from tkinter import ttk, IntVar
from tkinter.filedialog import askopenfilenames

# --- scientific Python ---
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# --- SciPy ---
import scipy
import scipy.fft as fft
from scipy import stats
from scipy.signal import find_peaks, savgol_filter, argrelextrema

# --- imaging ---
from skimage import io
from skimage.filters import threshold_otsu, threshold_local
from tifffile import imsave





# =========================
# Cross-correlation utility
# =========================
# Pearson cross-correlation over integer lags in [-maxShift, maxShift]
def xCorr(a, b, maxShift):

    corr = []
    for x in range (-maxShift, maxShift+1):  # iterate lag

        # circular shift (roll) and crop to suppress wrap-around artifacts
        cropA = np.roll(a, x)[maxShift:-maxShift]
        cropB = b[maxShift:-maxShift]

        # <<< added safeguard: check for empty arrays
        if cropA.size == 0 or cropB.size == 0:
            print(f"[xCorr] Warning: empty crop at shift {x}, maxShift={maxShift}, len(a)={len(a)}, len(b)={len(b)}")  # <<< changed
            corr.append(np.nan)
            continue

        if np.max(cropA) > 0 and np.max(cropB) > 0:
            corr.append(np.corrcoef(cropA, cropB)[0, 1])
        else:
            corr.append(np.nan)   #

    return(corr)


# ==================================================================
# Apply cross-correlation to CSV line profiles from NanoNet - Detect
# ==================================================================
# Args: path, file, chA/chB (0-based), pixelsize [nm], maxShift [px], chop flag/size [nm], randomizeShift [px]
# Returns: correlation matrix r(lag), piece counts, total lengths, chopped lengths
def doCorrelation(path, file, chA, chB, pixelsize, maxShift, chop, chopsizeNM, randomizeShift):# channels start at 0; returns correlations and piece counts
    global a
    global tableA

    allPieces = []
    allLength = []
    allLengthChop = []

    ch1 = pd.read_csv(path+file, index_col = 0)
    ch2 = pd.read_csv(path+file.replace("ch1", "ch2"), index_col = 0)
    ch3 = pd.read_csv(path+file.replace("ch1", "ch3"), index_col = 0)
    ch1 = ch1.set_index('Distance [nm]')
    ch2 = ch2.set_index('Distance [nm]')
    ch3 = ch3.set_index('Distance [nm]')

    dataset = [ch1, ch2, ch3]  # three channels

    tableA = dataset[chA]
    tableB = dataset[chB]

    correlations = []
    names = []
    for col in tableA.columns:  # iterate profiles (one per column)

        a = tableA[col].values
        if(randomizeShift>0):
            a = np.roll(a, randomizeShift)  # optional fixed circular lag

        b = tableB[col].values

        # restrict to overlapping support where signal > 0
        onsetA = np.where(a>0, np.linspace(0, len(a)+1, len(a), dtype = 'int'), np.linspace(0, len(a)+1, len(a), dtype = 'int')*0).min()
        endA = np.where(a>0, np.linspace(0, len(a)+1, len(a), dtype = 'int'), np.linspace(0, len(a)+1, len(a), dtype = 'int')*0).max()
        onsetB = np.where(b>0, np.linspace(0, len(b)+1, len(b), dtype = 'int'), np.linspace(0, len(b)+1, len(b), dtype = 'int')*0).min()
        endB = np.where(b>0, np.linspace(0, len(b)+1, len(b), dtype = 'int'), np.linspace(0, len(b)+1, len(b), dtype = 'int')*0).max()
        a = a[max(onsetA, onsetB):min(endA, endB)]# crop to overlapping nonzero range
        b = b[max(onsetA, onsetB):min(endA, endB)]# crop to overlapping nonzero range

        if(chop == True):

            chopsizePix = math.ceil(chopsizeNM/pixelsize) # segment length [px]
            pieces = math.floor(len(a)/chopsizePix)       # number of full segments

            count = 0
            for p in range(0, pieces):
                chopA = a[p*chopsizePix:(p+1)*chopsizePix]
                chopB = b[p*chopsizePix:(p+1)*chopsizePix]
                if(max(chopA)>0):
                    if(max(chopB)>0):
                        count+= 1

                if((len(chopA)>minSize) and (len(chopB)>minSize) ):
                    correlations+= [xCorr(chopA, chopB, maxShift)]  # r vs lag per segment
                    names+= [file+"_"+col+"_"+str(p)]

                allPieces+= [count]
            allLengthChop+= [count*chopsizeNM]# total analyzed length [nm] after chopping

        else:
            allPieces+= [1]
            allLengthChop+= [np.nan]

            if((len(a)>minSize) and (len(b)>minSize) ):
                correlations+= [xCorr(a, b, maxShift)]  # full-length correlation curve
                names+= [file+"_"+col+"_"]

        allLength+= [len(a)*pixelsize]# profile length [nm]

    correlations = pd.DataFrame(correlations).T  # rows: lags, cols: profiles/segments
    correlations.columns = names
    correlations = correlations.dropna(axis = 1) # remove empty/invalid series

    return correlations, allPieces, allLength, allLengthChop


# ============
# Peak finding
# ============
# Detect prominent local maxima per column; optional plotting
def findPeaks(inData, plot, height): # Peak detection per column; plot∈{True,False}; height = prominence threshold
    allMax = []
    columnHeaders = []

    for col in inData.columns:
        columnHeaders.append(col)

        tmpTrace = inData[col].values
        tmpTrace = scipy.ndimage.gaussian_filter(tmpTrace, 1.5, mode = 'reflect') # Gaussian smoothing (e.g. σ=1.5)

        maxima = find_peaks(tmpTrace, prominence = height, wlen = 15)[0].tolist() # indices of prominent local maxima
        maximaY = []
        for hit in maxima:
            maximaY+= [tmpTrace[hit]]  # peak amplitudes (not used later)

        maxima = np.array(maxima)
        maxima = np.ndarray.astype(maxima-(maxShift), int) # convert to lag domain by subtracting maxShift

        if(plot == True): # optional visualization
            fig = plt.figure()

        if (len(maxima)>0):
            allMax.append(maxima.tolist()) # store detected peak lags for this column

            if(plot == True):
                print(maxima)
                plt.scatter(maxima, tmpTrace[maxima], s = 600, c = 'red', marker = 'x')
        else:
            allMax.append([np.nan]) # no peaks detected

        if(plot == True):
            plt.plot(inData[col])                    # raw trace
            plt.plot(inData.index.values, tmpTrace)  # smoothed trace
            plt.title(col)
            plt.xlabel('pixels')
            plt.ylabel('intensity')
            plt.show()

    return allMax, columnHeaders


# =========================
# Histogram visualization(s)
# =========================
# Histogram of detected peak lags (nm); optional thresholding and first-peak selection
def plotHist(allMax, pixelsize, outDir, title, onlyFirst, removeBelowRes):
    fig = plt.figure()  # new figure
    maxBin = maxShift*pixelsize  # histogram upper bound [nm]
    global maxis
    maxis = pd.DataFrame(allMax)*pixelsize  # peaks → nm

    if(removeBelowRes == True):
        maxis = maxis[abs(maxis)>resolution]  # suppress sub-resolution lags

    if(onlyFirst == True):
        PeaksPx = abs(maxis).min(axis = 1)  # nearest-to-zero lag per sample
    else:
        PeaksPx = pd.DataFrame(np.ravel(maxis))  # flatten across columns

    PeaksPx.hist(bins = int(30/300*maxBin)*2, range = [0, maxBin], align = 'left')  # bin count ~ 0.2·maxBin
    plt.xlim(0, maxBin)
    #plt.title(title+"\n removeBelowRes = "+str(removeBelowRes))
    plt.xlabel("Distance [nm]")
    plt.ylabel("Counts")
    plt.savefig(outDir+"H_"+title+".png")  # write PNG

# Histogram of correlation coefficients (r ∈ [-1, 1]) for detected peaks
def plotHistStrength(allMax, pixelsize, outDir, title, onlyFirst, removeBelowRes):
    fig = plt.figure()                              # new figure
    maxBin = maxShift*pixelsize                     # computed but unused
    global maxis
    maxis = pd.DataFrame(allMax)                    # peak strengths (unitless)

    maxis.hist(bins = 25, range = [-1, 1], align = 'left')  # fixed r-domain binning
    plt.title(title+"\n removeBelowRes = "+str(removeBelowRes))  # annotate config
    plt.xlabel("Correlation")
    plt.ylabel("Counts")
    plt.savefig(outDir+"H_"+title+".png")           # write PNG

    plt.savefig(outDir+"C_"+title+".png")  # write PNG

# Correlation trace
def plotCorr(matrix, pixelsize, outDir, title):
    fig = plt.figure()  # new figure
    plt.plot(correlations.index.values*pixelsize, np.mean(correlations, axis = 1))  # mean r across columns vs physical lag (nm)
    plt.xlabel("Distance [nm]")  # lag in nanometers
    plt.ylabel("Correlation")   # Pearson r
    #plt.title(title)            # dataset label
    plt.savefig(outDir+"C_"+title+".png")  # write PNG
    plt.close('all')


# =============
# Parameter GUI
# =============
global chop           # segment profiles if True
global chopsizeNM     # window length [nm]
global pixelsize      # sampling [nm/px]
global maxShift       # maximal shift [px]
global plotRawData_min  # debug: anti-correlation traces
global plotRawData_max  # debug: correlation traces
global onlyFirst      # keep only first detected peak if True

global path
global files

root = Tk()  # Tk root for file selection

PathList = askopenfilenames(filetypes = (
(".csv", "*_ch1.csv"), ("All files", "*.*")
))

root.destroy()  # close file dialog/root

files = []
for thisfile in PathList:
    folder, file = os.path.split(thisfile)
    files = files+[file]

path = folder+"//"  # working directory

# defaults (overwritable via GUI)
chop = True  # segment into 'chopsizeNM' windows
chopsizeNM = 1500  # [nm]

try:
    pixelsize = round(float(pd.read_csv(os.path.join(path, files[0]), header=None).iloc[2, 1])*1000)
except Exception as e:
    print(f"Warning: could not read pixelsize from {files[0]} (B3). Using default = 1. Error: {e}")
    pixelsize = 1

maxShift = 40      # [px]
plotRawData_min = False  # debug plot: anti-correlation
plotRawData_max = False  # debug plot: correlation
onlyFirst = True         # only first (anti-)correlation peak
minSize = 30             # minimal profile length [px]
randomizeShift = 0       # optional fixed circular shift [px]
height = 0.03             # peak prominence threshold (0–2)
resolution = 50        # microscope resolution [nm]
ignoreZeroPeak = True   # False/True: ignore sub-resolution peaks in AUTO-correlation histograms

def getParameters():
    # Simple GUI to review/override parameters; writes back to globals
    def show_entry_fields():
        global chop, chopsizeNM, pixelsize, maxShift, plotRawData_min, plotRawData_max
        global onlyFirst, minSize, randomizeShift, height, resolution, ignoreZeroPeak

        chop             = v.get()
        chopsizeNM       = e2.get()
        pixelsize        = e3.get()
        maxShift         = e4.get()
        plotRawData_min  = v5.get()
        plotRawData_max  = v6.get()
        onlyFirst        = v7.get()
        minSize          = e8.get()
        randomizeShift   = e9.get()
        height           = e10.get()
        resolution       = e11.get()
        ignoreZeroPeak   = v12.get()

        master.destroy()  # close parameter window

    def gogogo():
        global addChannel
        addChannel = False
        show_entry_fields()

    master = Tk()
    master.title("NanoNET:Analysis")

    # ---- labels ----
    Label(master, text="                                       ").grid(row=0,  column=0, sticky=W)
    Label(master, text="Break profile into fragments           ").grid(row=1,  column=0, sticky=W)
    Label(master, text="Fragment size [nm]:             ").grid(row=2,  column=0, sticky=W)
    Label(master, text="                                       ").grid(row=3,  column=0, sticky=W)
    Label(master, text="Pixelsize [nm]:                             ").grid(row=4,  column=0, sticky=W)
    Label(master, text="Shift Range [pixels]:               ").grid(row=5,  column=0, sticky=W)
    Label(master, text="                                       ").grid(row=6,  column=0, sticky=W)
    Label(master, text="Debugging: Plot raw anticorrelation data          ").grid(row=7,  column=0, sticky=W)
    Label(master, text="Debugging: Plot raw correlation data              ").grid(row=8,  column=0, sticky=W)
    Label(master, text="                                       ").grid(row=9,  column=0, sticky=W)
    Label(master, text="Plot only first (anti-)correlation peak").grid(row=10, column=0, sticky=W)
    Label(master, text="                                       ").grid(row=11, column=0, sticky=W)
    Label(master, text="Minimal profile length [pixels]        ").grid(row=12, column=0, sticky=W)
    Label(master, text="                                       ").grid(row=13, column=0, sticky=W)
    Label(master, text="Randomize by shift [pixels]            ").grid(row=14, column=0, sticky=W)
    Label(master, text="                                       ").grid(row=15, column=0, sticky=W)
    Label(master, text="Minimal height for peak finding (0-2)  ").grid(row=16, column=0, sticky=W)
    Label(master, text="                                       ").grid(row=17, column=0, sticky=W)
    Label(master, text="Microscope resolution [nm]").grid(row=18, column=0, sticky=W)
    Label(master, text="Suppress sub-resolution peaks in AUTO-corr").grid(row=19, column=0, sticky=W)

    Label(master, text="                                       ").grid(row=20, column=0, sticky=W)
    Label(master, text="Please cite Winkler et al. 2025. Methodical details: https://github.com/ngimber/NanoNET").grid(row=21, column=0, sticky=W)

    # ---- inputs ----
    v = IntVar(value=chop)
    Checkbutton(master, variable=v).grid(row=1, column=2, sticky=W)

    v2 = IntVar(value=chopsizeNM)
    e2 = Entry(master, text=v2, width=45); e2.grid(row=2, column=2, sticky=W)

    v3 = IntVar(value=pixelsize)
    e3 = Entry(master, text=v3, width=45); e3.grid(row=4, column=2, sticky=W)

    v4 = IntVar(value=maxShift)
    e4 = Entry(master, text=v4, width=45); e4.grid(row=5, column=2, sticky=W)

    v5 = IntVar(value=plotRawData_min)
    Checkbutton(master, variable=v5).grid(row=7, column=2, sticky=W)

    v6 = IntVar(value=plotRawData_max)
    Checkbutton(master, variable=v6).grid(row=8, column=2, sticky=W)

    v7 = IntVar(value=onlyFirst)
    Checkbutton(master, variable=v7).grid(row=10, column=2, sticky=W)

    v8 = IntVar(value=minSize)
    e8 = Entry(master, text=v8, width=45); e8.grid(row=12, column=2, sticky=W)

    v9 = IntVar(value=randomizeShift)
    e9 = Entry(master, text=v9, width=45); e9.grid(row=14, column=2, sticky=W)

    v10 = IntVar(value=height)
    e10 = Entry(master, text=v10, width=45); e10.grid(row=16, column=2, sticky=W)

    # NEW fields:
    v11 = IntVar(value=resolution)
    e11 = Entry(master, text=v11, width=45); e11.grid(row=18, column=2, sticky=W)

    v12 = IntVar(value=ignoreZeroPeak)
    Checkbutton(master, variable=v12).grid(row=19, column=2, sticky=W)

    Button(master, text='Start Processing', fg="green", bg="gray83",
           command=gogogo).grid(row=21, column=4, columnspan=1000)

    master.attributes("-topmost", True)
    mainloop()


getParameters()  # launch parameter GUI

# cast GUI-returned values to proper types
chopsizeNM = int(chopsizeNM)
pixelsize = int(pixelsize)
maxShift = int(maxShift)
minSize = int(minSize)
randomizeShift = int(randomizeShift)
height = float(height)
height = float(height)
resolution = float(resolution)
ignoreZeroPeak = bool(int(ignoreZeroPeak))

# echo configuration to console
print(files)
print(folder)
print(chop)
print(chopsizeNM)
print(pixelsize)
print(maxShift)
print(plotRawData_min)
print(plotRawData_max)
print(onlyFirst)
print(minSize)
print(randomizeShift)
print(height)



# ====================
# Output dir structure
# ====================

# Define base directories
plotDir = os.path.join(path+"results//", "plots")
summaryplotDir = os.path.join(path+"results//", "summaryplots")
tableDir = os.path.join(path+"results//", "tables")
logDir = os.path.join(path+"results//", "log_files")

# Define subfolders
plot_subdirs = ["histograms", "correlograms"]
table_subdirs = ["correlograms", "maxima","stats"]

def ensure_subdirs(base_dir, subdirs):
    # always make the base dir
    os.makedirs(base_dir, exist_ok=True)
    print(f"Ensured base directory: {base_dir}")

    # then handle subdirs if any
    for sub in subdirs:
        folder = os.path.join(base_dir, sub)
        if os.path.exists(folder):
            print(f"Exists: {folder}")
        else:
            os.makedirs(folder)
            print(f"Created: {folder}")


# Make plot-related folders
ensure_subdirs(plotDir, plot_subdirs)
ensure_subdirs(summaryplotDir, plot_subdirs)
ensure_subdirs(tableDir, table_subdirs)
ensure_subdirs(logDir,[])


# ==========
# Write Log Files
# ==========
pd.DataFrame([[files], [folder], [chop], [chopsizeNM], [pixelsize], [maxShift], [onlyFirst], [minSize], [randomizeShift],[resolution],[ignoreZeroPeak]], index = ["files", "folder", "chop", "chopsizeNM", "pixelsize", "maxShift", "onlyFirst", "minSize", "randomizeShift","resolution","ignore zero peak"]).to_csv(logDir+"//log"+datetime.now().strftime("%Y%m%d_%H%M%S")+".txt")


# ==========
# Processing
# ==========
import matplotlib.pyplot as plt
import pandas as pd

# set plotting style
plt.style.use("seaborn-v0_8-dark")

combis = []

for chA in range(3):
    for chB in range(3):
        print(f"channel A: {chA}", f"channel B: {chB}")
        c = f"{min(chA, chB)}_{max(chA, chB)}"

        if c in combis:
            print(f"Skipping duplicate combination: {c}")
            continue
        else:
            combis.append(c)

        allallMax = []
        nChops = []
        profileLength = []
        peaksPerMicron = []
        profileLengthChop = []
        result_list = []

        for f, file in enumerate(files):
            title = f"{file}_{chA}_{chB}"

            # compute correlation
            corr = doCorrelation(path, file, chA, chB, pixelsize, maxShift, chop, chopsizeNM, randomizeShift)
            correlations = corr[0]
            nChops.append(corr[1])
            profileLength.append(corr[2])            # [nm]
            profileLengthChop.append(corr[3])        # [nm]
            correlations.index = correlations.index - maxShift

            # detect peaks
            allMax, columnHeaders = findPeaks(correlations, plotRawData_max, height)
            peak_density = pd.DataFrame(allMax).count(axis=1).values / (2 * maxShift * pixelsize / 1000)
            peaksPerMicron.append(peak_density.tolist())


            allallMax += allMax
            correlations.to_csv(f"{tableDir}/correlograms/corr_{file[:-8]}_{chA}_{chB}.csv")
            pd.DataFrame(allMax).T.to_csv(f"{tableDir}/maxima/summary_detected_maxima_pix_{file[:-8]}_ch{chA}_ch{chB}.csv")
            pd.DataFrame(allallMax).to_excel(f"{tableDir}/maxima/detected_maxima_pix_ch{chA}_ch{chB}.xlsx")



            if(bool(ignoreZeroPeak)==True):
                removeBelowRes = (chA == chB)#remove peaks below reslution on (only for auto-correlation)
            else:
                removeBelowRes=False

            # per-file plots
            plotHist(allMax, pixelsize, plotDir+"//histograms//", f"MAX{file[:-4]}_{chA}_{chB}", onlyFirst, removeBelowRes)
            plotCorr(correlations, pixelsize, plotDir+"//correlograms//", f"{file[:-4]}_{chA}_{chB}")

            # collect correlation data
            if f == 0:
                allCorrelations = correlations
            else:
                allCorrelations = pd.concat([allCorrelations, correlations], axis=1)

            result_list.append(allCorrelations)

        # final averaging and summaries
        result_df = pd.concat(result_list, axis=1)
        av_CorrRow = result_df.mean(axis=1)

        # Save correlogram summaries
        #result_df.to_csv(f"{tableDir}/correlograms/corr_summary_ch{chA}_ch{chB}.csv")

        # Save correlation averages
        pd.DataFrame(av_CorrRow).to_excel(f"{tableDir}/correlograms/corr_average_ch{chA}_ch{chB}.xlsx")

        # Save peaks per micron
        pd.DataFrame(peaksPerMicron, index=files).T.to_csv(
            f"{tableDir}/stats/corr_PeaksPerMicronShift_ch{chA}_ch{chB}.csv"
)


        # summary plots
        plotCorr(allCorrelations, pixelsize, summaryplotDir+"//correlograms//", f"summary_ch{chA}_ch{chB}")
        plotHist(allallMax, pixelsize, summaryplotDir+"//histograms//", f"summary_MAX_ch{chA}_ch{chB}", onlyFirst, removeBelowRes)

    # per-channel metrics
    pd.DataFrame(nChops, index=files).T.to_excel(f"{tableDir}/stats/ChoppedPieces.xlsx")
    pd.DataFrame(profileLength, index=files).T.to_excel(f"{tableDir}/stats/ProfileLengthInNM.xlsx")
    #pd.DataFrame(profileLengthChop, index=files).T.to_excel(f"{tableDir}/stats/ProfileLengthChopInNM.xlsx")
