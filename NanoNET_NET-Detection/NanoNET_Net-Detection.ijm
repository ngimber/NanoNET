// ========================================================
//  User Parameters
// ========================================================

linewidthUM          = 0.21;// Width of skeleton lines and profile traces (µm)
minlengthUM          = 0.225;// Minimum length for any measured profile (µm)
dnaChannel         = 3;// DNA image channel to process (1-based index)
dir = getDirectory("home");
extention = ".tif";

// ========================================================
//  Contrast Settings
// ========================================================

// Flag: correct local contrast if true
locContrast = true;

// Junction handling method: "NONE" for splits, "SLOPE" for connections
IntersectMethod = "NONE";


// ========================================================
//  Batch Mode & Unified Dialog Setup
// ========================================================

// Enable batch processing
setBatchMode(true);

// Prepare dialog for all inputs (including folder selection)
Dialog.create("NanoNET:Detection");

// --- File & Folder Settings ---
Dialog.addMessage("***  File & Folder Settings  ***");
// Select analysis folder
Dialog.addDirectory("Data folder:", dir);
Dialog.addMessage("\n");

// --- Import Parameters ---
Dialog.addMessage("***  Import Parameters  ***");
// File extension filter
Dialog.addString("File extension:", extention);
Dialog.addMessage("\n");

// --- Filter Profiles ---
Dialog.addMessage("***  Filter Profiles  ***");
Dialog.addNumber("Minimal profile length (µm):", minlengthUM);
Dialog.addMessage("\n");

// --- Segmentation Parameters ---
Dialog.addMessage("***  Segmentation Parameters  ***");
Dialog.addNumber("DNA channel (first = 1):", dnaChannel);
Dialog.addCheckbox("Correct local contrast:", locContrast);
Dialog.addMessage("\n");

// --- Measurement Parameters ---
Dialog.addMessage("***  Measurement Parameters  ***");
Dialog.addNumber("Line width (µm):", linewidthUM);
Dialog.addMessage("\n");

Dialog.show();

// Fetch dialog inputs
dir           = Dialog.getString();
if(endsWith(dir, "/")==false){dir=dir+"\\";}

extention   = Dialog.getString();
minlengthUM   = Dialog.getNumber();
dnaChannel  = Dialog.getNumber();
locContrast = Dialog.getCheckbox();
linewidthUM   = Dialog.getNumber();

// Retrieve list of files in folder
files = getFileList(dir);


// ========================================================
//  Logging Setup
// ========================================================

print("Minimal profile length: " + minlengthUM);
print("DNA channel: " + dnaChannel);
print("Correct local contrast: " + locContrast);
print("Line width (µm): " + linewidthUM);
selectWindow("Log");

// Save log file with timestamp
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
saveAs("txt", dir + "log_" + year + month + dayOfMonth + "_" + hour + minute + ".txt");
print("\\Clear");




// ========================================================
//  Main Processing Loop
// ========================================================
for (i = 0; i < files.length; i++) {
    if (indexOf(files[i], extention) >= 0) {

        // ----------------------------------------------------
        //  Cleanup: close images and free memory
        // ----------------------------------------------------
        while (nImages > 0) { selectImage(nImages); close(); }
        run("Collect Garbage"); call("java.lang.System.gc");
     

        // ----------------------------------------------------
        //  Open image and wait for ready
        // ----------------------------------------------------
        open(dir + files[i]);
        while (getTitle != files[i]) { wait(500); print("wait to load files"); }
        getPixelSize(unit, pixelWidth, pixelHeight);
        minlength   = Math.ceil(minlengthUM/pixelWidth); //convert into pixels
		linewidth   = Math.ceil(linewidthUM/pixelWidth); //convert into pixels
		

        // ----------------------------------------------------
        //  Filter Parameters
        // ----------------------------------------------------
        
		
		// Sigma for Non-local Means Denoising filter
		denoiseSigma       = Math.ceil(0.45/pixelWidth); 
		
		// Radius for disk-shaped Closing filter 
		diskClosingRadius   = Math.ceil(0.075/pixelWidth); 
		
		// Radius for square-shaped Closing filter when local contrast is off 
		squareClosingRadius = Math.ceil(0.075/pixelWidth); 
		
		// Minimum particle size before skeletonization 
		preSkeletonSize    = Math.ceil(0.15/pixelWidth); 
		

        // ----------------------------------------------------
        //  Setup output directories
        // ----------------------------------------------------
        title = getTitle();
        dir   = getDirectory("image");
        subfolder = dir + "\\analysis\\";
        File.makeDirectory(subfolder);
        File.makeDirectory(subfolder + "ROIs\\");
        File.makeDirectory(subfolder + "CSVs\\");
        File.makeDirectory(subfolder + "CSVs_param\\");
        File.makeDirectory(subfolder + "Skeleton\\");

        // ----------------------------------------------------
        //  ROI Manager reset & duplicate DNA channel
        // ----------------------------------------------------
        roiManager("reset");
        run("Select All");
        run("Duplicate...", "title=tmp duplicate channels=" + dnaChannel);

        // ----------------------------------------------------
        //  Local contrast correction (optional)
        // ----------------------------------------------------
        if (locContrast) {
            run("Normalize Local Contrast", "block_radius_x=1 block_radius_y=1 standard_deviations=30 center stretch");
            setAutoThreshold("Otsu dark");
            run("Convert to Mask");
        }

        // ----------------------------------------------------
        //  Denoising & morphological closing
        // ----------------------------------------------------
        run("Non-local Means Denoising", "sigma=" + denoiseSigma + " smoothing_factor=1");
        run("Morphological Filters", "operation=Closing element=Disk radius=" + diskClosingRadius);

        // ----------------------------------------------------
        //  Ridge detection if no local contrast
        // ----------------------------------------------------
        if (!locContrast) {
            run("Top Hat...", "radius=1");
            run("Select All"); run("Apply LUT"); run("8-bit");
            run("Ridge Detection", "line_width=1 high_contrast=230 low_contrast=0.1 extend_line displayresults add_to_manager make_binary method_for_overlap_resolution=" + IntersectMethod + " sigma=1.1 lower_threshold=0.1 upper_threshold=8 minimum_line_length=0 maximum=0");
            run("Morphological Filters", "operation=Closing element=Square radius=" + squareClosingRadius);
        }

        // ----------------------------------------------------
        //  Skeletonize and skeleton analysis
        // ----------------------------------------------------
        roiManager("reset");
        run("Skeletonize (2D/3D)");
        run("Analyze Skeleton (2D/3D)", "prune=none calculate show display");
        wait(500);
        selectWindow("Tagged skeleton");
        setThreshold(126, 128);
        run("Analyze Particles...", "size=" + preSkeletonSize + "-Infinity pixel clear add");
        roiManager("Combine"); run("Clear Outside"); run("Convert to Mask");

        // Set line settings and detect final skeleton particles
        setLineWidth(linewidth);
        setAutoThreshold("Default dark no-reset");
        run("Analyze Particles...", "size=" + minlength + "-Infinity pixel clear add");

        // Close temp windows and save skeleton image
        close("tmp-Closing-Closing-labeled-skeletons");
        close("Longest shortest paths");
        close("tmp-Closing-Closing");
        close("tmp-Closing");
        close("tmp-Closing Detected segments");
        close("tmp");
        saveAs("Tiff", subfolder + "Skeleton\\Skeleton_" + title);

        // ========================================================
        //  Profile Measurement
        // ========================================================
        run("Plots...", "width=450 height=200 font=14 auto-close list minimum=0 maximum=0");
        selectImage(title);
        getDimensions(width, height, channels, slices, frames);
        getPixelSize(unit, pixelWidth, pixelHeight);
        roiManager("Save", subfolder + "ROIs\\raw_" + title + ".zip");
        nRois = roiManager("count");

        for (ch = 1; ch <= channels; ch++) {
            // Clear previous results and initialize distance column
            run("Clear Results"); setResult("Distance [nm]", 0, 0);
            selectImage(title); setSlice(ch);
		
		setBatchMode(false);
		for (n=0;n<nRois;n++)
			{
			roiManager("Select", n);
			Roi.setStrokeWidth(linewidth);
			run("Area to Line");
			roiManager("Add");
			//run("Plot Profile");
			profile=getProfile();
			writeResulttabe(profile, n);			
			}
		setBatchMode(true);
            // Convert pixel indices to nanometer distances
            for (row = 0; row < nResults; row++) {
                setResult("Distance [nm]", row, row * pixelWidth);
            }
            close("Plot Values");
            saveAs("Results", subfolder + "CSVs\\Profile_" + title + "_ch" + ch + ".csv");
        }

        // ========================================================
        //  ROI Cleanup & Additional Measurements
        // ========================================================
        // Save finalized line ROIs
        for (n = 0; n < nRois; n++) { roiManager("Select", 0); roiManager("Delete"); }
        roiManager("Save", subfolder + "ROIs\\line_" + title + ".zip");

        // Set detailed measurement options
        run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack display redirect=None decimal=9");
        run("Clear Results"); selectImage(title);
        roiManager("Deselect"); roiManager("multi-measure measure_all");
        saveAs("Results", subfolder + "CSVs_param\\Profile_" + title + ".csv");

        // ========================================================
        //  Line Thickness Measurement
        // ========================================================
        selectImage(title);
        roiManager("reset"); roiManager("Open", subfolder + "ROIs\\raw_" + title + ".zip");
        setLineWidth(1);
        nRois = roiManager("count");
        for (n = 0; n < nRois; n++) { roiManager("Select", n); run("Area to Line"); roiManager("Add"); }
        for (n = 0; n < nRois; n++) { roiManager("Select", 0); roiManager("Delete"); }
        run("Clear Results");
        run("Set Measurements...", "mean stack display redirect=None decimal=9");
        run("Duplicate...", "title=thr duplicate channels=3");
        setOption("BlackBackground", true);
        setAutoThreshold("MinError"); run("Convert to Mask"); run("Convert to Mask");
        run("Distance Map");
        roiManager("Deselect"); roiManager("multi-measure measure_all");
        saveAs("Results", subfolder + "CSVs_param\\Linewidth_" + title + ".csv");

        // ----------------------------------------------------
        //  Final Cleanup: close summary windows
        // ----------------------------------------------------
        close("Branch information"); close("Junctions"); close("*");
    }
}

// ========================================================
//  Helper Functions
// ========================================================

// Reads a result column into a Java array
function readArray(columnname) {
    storageArray = newArray(nResults);
    for (row = 0; row < nResults; row++) {
        storageArray[row] = getResult(columnname, row);
    }
    return storageArray;
}

// Writes an array into the Results table under given column name
function writeResulttabe(input, columnname) {
    len = lengthOf(input);
    for (row = 0; row < len; row++) {
        setResult(columnname, row, input[row]);
    }
}

exit;
