# NanoNET Toolbox  
**Generation and Analysis of NET Filament Line Profiles**

This repository provides tools for studying neutrophil extracellular traps (NETs) from light microscopy data.  
NanoNET consists of the following workflows:

- **NanoNET Detection** – Automated segmentation and extraction of NET filament profiles from light microscopy data (e.g., SIM, STED, SMLM).  
- **NanoNET Analysis** – Statistical correlation and periodicity analysis of filament profiles.

Detailed information on sample preparation, imaging, and NanoNET can be found in the bioRxiv paper:  
**Nanoscale Mapping Reveals Periodic Organization of Neutrophil Extracellular Trap Proteins** — bioRxiv (2025).  
DOI: https://doi.org/10.1101/2025.07.28.665103

---

## 📦 Workflows

### 1. NET Detection
**Purpose:** Extract NET filament line profiles from 2D/3D microscopy stacks. Line profiles for all channels (three-channel images expected) are measured along a continuous DNA backbone channel.

![NanoNET Detection GUI](https://github.com/ngimber/NanoNET/blob/main/NanoNET_NET-Detection/GUI_NanoNET-Detection.png)

**Workflow steps:**
1. **Preprocessing**
   - (3D only) Create **maximum-intensity projections** from z-ranges containing NET signal.
   - (Optional) **Cropping:** remove non-NET regions or artifacts (e.g., cell remnants).

2. **Automated segmentation (NanoNET Detection)**  
   **Set parameters:**
   - **Data folder:** select the folder containing TIFF files (2D images or 3D projections).
   - **File extension:** choose the image file type (e.g. .tif).
   - **Minimal profile length (µm):** profiles shorter than this threshold are ignored.
   - **Channel selection:** define the **DNA backbone** channel (requires continuous DNA labeling).
   - **Correct local contrast:** compensates for shading/illumination inhomogeneity.
   - **Line width (px):** width over which intensities are integrated for line profiles.

**Outputs**
- Skeletonized NET filaments  
- Per-channel line intensity profiles (CSV)  
- Parameter/configuration files for reproducibility

---

### 2. NET Analysis
**Purpose:** Perform correlation analysis of NET filament line profiles to identify structural periodicities and colocalization patterns.

![NanoNET Analysis GUI](https://github.com/ngimber/NanoNET/blob/main/NanoNET_NET-Analysis/GUI_NanoNET-Analysis.png)

**Workflow steps:**
1. **Input:** line intensity profiles from **NET Detection** (CSV files).

2. **Automated analysis**  
   **Set parameters:**
   - **Fragments:** enable **Break into fragments** to split profiles into segments of defined length (e.g., 1500 nm).
   - **Pixel size (nm/px):** set the image sampling.
   - **Resolution:** expected resolution of the microscopy technique.
   - **Plot only first peak:** if enabled, only the first peak in the correlogram is considered.
   - **Shift intervals:** maximum lags for auto-/cross-correlation (steps are pixels).
   - **Minimal filament length:** exclude short fragments from analysis.

**Outputs**
- Correlation tables (auto- & cross-correlation values)  
- Averaged correlation profiles  
- Extracted periodicities (histograms)

---

## ⚙️ Requirements

- **NET Detection**
  - Fiji/ImageJ with NanoNET plugin — https://imagej.net/software/fiji/  

- **NET Analysis**
  - Python 3.9.13  
  - Jupyter Notebook
### 📥 Python package installation

Install the Python dependencies with `pip` (recommended in a virtual environment):

```bash

# Upgrade pip
pip install --upgrade pip

# Core scientific stack
pip install pandas numpy scipy matplotlib

# Image IO and processing
pip install scikit-image tifffile imageio

# Jupyter (for running the notebooks/GUI)
pip install notebook ipykernel

# Optional: plotting helpers
pip install seaborn
```

---

## 📖 References

1. **NanoNET (study)**  
   Gimber N., *et al.* (2025). **Nanoscale Mapping Reveals Periodic Organization of Neutrophil Extracellular Trap Proteins.** *bioRxiv*. https://doi.org/10.1101/2025.07.28.665103

2. **Image processing methods**  
   Otsu, N. (1979). **A threshold selection method from gray-level histograms.** *IEEE Transactions on Systems, Man, and Cybernetics*, 9(1), 62–66. https://doi.org/10.1109/TSMC.1979.4310076  
   Buades, A., Coll, B., & Morel, J.-M. (2005). **A non-local algorithm for image denoising.** *CVPR 2005*, 2, 60–65. https://doi.org/10.1109/CVPR.2005.38

3. **ImageJ / Fiji**  
   Schindelin, J., *et al.* (2012). **Fiji: an open-source platform for biological-image analysis.** *Nature Methods*, 9(7), 676–682. https://doi.org/10.1038/nmeth.2019  
   Schneider, C. A., Rasband, W. S., & Eliceiri, K. W. (2012). **NIH Image to ImageJ: 25 years of image analysis.** *Nature Methods*, 9(7), 671–675. https://doi.org/10.1038/nmeth.2089

4. **Python ecosystem**  
   Virtanen, P., *et al.* (2020). **SciPy 1.0: fundamental algorithms for scientific computing in Python.** *Nature Methods*, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2  
   Harris, C. R., *et al.* (2020). **Array programming with NumPy.** *Nature*, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2  
   Hunter, J. D. (2007). **Matplotlib: A 2D graphics environment.** *Computing in Science & Engineering*, 9(3), 90–95. https://doi.org/10.1109/MCSE.2007.55  
   McKinney, W. (2010). **Data Structures for Statistical Computing in Python.** *Proc. of the 9th Python in Science Conf. (SciPy 2010)*, 51–56. https://doi.org/10.25080/Majora-92bf1922-00a  
   van der Walt, S., *et al.* (2014). **scikit-image: image processing in Python.** *PeerJ*, 2, e453. https://doi.org/10.7717/peerj.453  
   Gohlke, C. (tifffile). **tifffile: Read and write TIFF files.** https://pypi.org/project/tifffile/
