//This macro opens a whole folder of Incucyte data one by one (looks for hyphen, not underscore to split parts of filename for parsing and looping), turns each into a time series, 
///and calculates the earliest timepoint any GFP+ cells were found in any colonies detected in the final frame larger than 20pixels
/// and saves ROIs, results, the stack of masks for the full time series, and closes everything
//9/9/26 Anneliese Jarman, Imaging Facility, Babraham Institute

// Select input folder and get list of files
inputFolder = getDirectory("Please select the folder that contains your images"); 
outputFolder = getDirectory("Please select the folder you want your outputs to go into.");
fileList = getFileList(inputFolder);                                                         

// Create output folder
File.makeDirectory(inputFolder + "tmp");                                                 
tmp = inputFolder + File.separator + "tmp";
out1 = File.getParent(tmp);
File.delete(tmp); 

// Count GFP files so we can create the right sized array
nGFP = 0;
for (i = 0; i < fileList.length; i++) {
    if (matches(fileList[i], ".*_GFP_.*"))
        nGFP++;
}


// Store GFP filenames in an array
GFPfiles = newArray(nGFP);
index = 0;
for (i = 0; i < fileList.length; i++) {
    if (matches(fileList[i], ".*_GFP_.*")) {
        GFPfiles[index] = fileList[i];
        index++;
    }
}


// Sort filenames so timepoints are in chronological order
Array.sort(GFPfiles);


// Loop through all GFP files and identify each new prefix/well/FOV
filename = "";
for (i = 0; i < GFPfiles.length; i++) {

    // Extract prefix, channel, well, FOV and timepoint
    parts = split(GFPfiles[i], "-");

    // Ignore the timepoint when constructing the sequence name
    currentFilename = parts[0] + "-" + parts[1] + "-" ;


    // Only open a new sequence when the prefix/well/FOV changes
    if (currentFilename != filename) {

        filename = currentFilename;

        print("Opening sequence: " + filename);

        // Open all timepoints belonging to this prefix/well/FOV as a time series
        File.openSequence(inputFolder, "filter=" + filename + " sort");
        rename("timeStack");

        ////PROCESSING and MASK MAKING////
        
        //duplicate to use one for mask and reserve one for redirected measurements (called filename)
        run("Duplicate...", "duplicate");
        // blur to smooth, and autothreshold to generate mask
		run("Gaussian Blur...", "sigma=2 stack");
		setOption("ScaleConversions", true);
		run("Auto Threshold", "method=Triangle white stack");
		rename("maskStack");
		
		////TRIMMING THE FRAMES BEFORE "BLOWOUTS" -if no foreground can be detected (no foreground there yet or not enough contrast to detect it with the autothreshold)////

		selectWindow("maskStack");
						
		//// Calculate foreground area for all slices in stack////
		f = nSlices; //frames 
		
		foregroundAreaStack = newArray(f+1);
		
		for (j = 1; j <= f; j++) {
		    setSlice(j);
		    getHistogram(values, counts, 256);
		    foregroundArea = counts[255];
		    foregroundAreaStack[j] = foregroundArea; //contains the foreground area for every frame in an array
		}
		
		keepIndex = 1; //if there are no blowOuts 
		for (k = f; k>=2; k--){
			if (foregroundAreaStack[k-1] >= 2*foregroundAreaStack[k]){
			//Cut off the bad slice from the start of the maskStack,  they have no foreground)
			keepIndex = k;
			break
			}
		}
		
		validFrames = f - keepIndex + 1 ;
		removedFrames = f - validFrames;
		
		print("f: " + f);
		print("keepIndex: " + keepIndex);
		print("validFrames: " + validFrames);
		print("removedFrames:" + removedFrames);
		print("substack: " + keepIndex + "-" + f);
		

		if ( validFrames <= 3) {
			print("Not enough valid frames. Series skipped:" + filename);
				selectWindow("Log");
				saveAs("Text", outputFolder + filename + "_SKIPPED_Log.txt");
				close("Log");
				close("*");
				continue;
		}	
		
		run("Make Substack...", "slices=" + keepIndex + "-" + f);
		rename("useableMaskStack");
		//Do the same for the data stack
		selectImage("timeStack");
		run("Make Substack...", "slices=" + keepIndex + "-" + f);
		rename("useableDataStack");
		
		//// Find the earliest timepoint of GFP+ cells within the GFP+ colonies of the final frame ////
		new_f = nSlices; 		
		//		selectWindow("useableMaskStack");
		selectWindow("useableMaskStack");
		run("Duplicate...", "duplicate range=" + new_f);
		rename("finalTimepointMask");
	
		//// Walk backwards through the REMAINING mask stack and encode frame numbers ====
		selectImage("useableMaskStack");            // stack has shrunk -- this is the new (smaller) slice count
		placeholder = f + 1;      // still based on the ORIGINAL f, so it can't collide with any real frame number
		
		for (l = new_f; l >= 1; l--) {
		    setSlice(l);
		    changeValues(0, 0, placeholder);   // background -> "no data" (change values between 0 and 0 to a placeholder number)
		    changeValues(255, 255, l);         // foreground -> this frame's number (change values between 255 and 255 to its frame number)
		}
		
		// Each slice is now labelled with either its own frame number (foreground) or the placeholder value (background)
		rename("encodedUseableMaskStack");
		
		// Next step: Z-project this time stack using MINIMUM intensity. // For each pixel, the minimum across all slices = the earliest frame number it was ever foreground in (or placeholder (number of frames + 1), if it was never reliably foreground)
		selectImage("encodedUseableMaskStack");
		run("Z Project...", "projection=[Min Intensity]");
		rename("earliestTimepointImg"); //this image encodes the earliest frame in which each pixel was classified as foreground
		
		// Within the final timepoint colony areas, find the minimum pixel value, encoding the first frame in which any pixels within this colony were detected! :D
		run("Set Measurements...", "min redirect=earliestTimepointImg decimal=3");
		selectWindow("finalTimepointMask");
		run("Analyze Particles...", "size=20-Infinity add display");
		
		//update the results table to get rid of the max measurement that we don't want, and rescale the min value to the original number of frames
		for (row = 0; row < nResults; row++) {
   			 minVal = getResult("Min", row);
    		setResult("EarliestFrame", row, minVal + removedFrames);
		}
		updateResults();	
		Table.deleteColumn("Max");
		Table.deleteColumn("Min");

		
		////SAVE AND CLOSE STUFF READY FOR THE NEXT ROUND////
	    
	    //save ROIs - Their order agrees with the results rows
		roiManager("save", outputFolder + filename + "_MeasuredColonies.zip");
		roiManager("reset");
		close("ROI Manager");
		
		//save images/masks of interest
		selectImage("maskStack");
		saveAs("Tiff", outputFolder + filename + "_maskStack");
		
		selectImage("earliestTimepointImg");
		saveAs("Tiff", outputFolder + filename + "_earliestTimepointImg");		

		//save and close results
		saveAs("Results", outputFolder + filename + "_Results.csv");
		close("Results");
		
		//save log. This tells you which the filenumber of the folder is being processed, its name, frames, the first frame we're keeping 
		//from (post-blowout), the total number of valid remaining frames, and which they are
		selectWindow("Log");
		saveAs("Text", outputFolder + filename + "_Log.txt");
		close("Log");
		
		//close images
		close("*");
		
    }
}

