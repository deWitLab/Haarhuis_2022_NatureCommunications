#@ File (label = "Input file", style = "file") input_file
#@ File (label = "Output folder", style = "directory") output

#@ Integer (label = "Nuclei channel", value=2, style="spinner", min=1, max=2) nucleiChannel
#@ Integer (label = "Signal of interest channel", value=1, style="spinner", min=1, max=2) signalChannel

#@ Integer (label = "Median filter radius before nuclei segmentation", value=1, style="spinner", min=0, max=10, value=4) medianFilterRadius
#@ Boolean (label = "Automatic nuclei thresholding?") autoThresholdNuclei
#@ Integer (label = "Lower nucleus diameter limit (um)", style = "spinner", min=0, max=1000, value=4) lowerDiameterLimit
#@ Integer (label = "Upper nucleus diameter limit (um)", style = "spinner", min=0, max=1000, value=40) upperDiameterLimit
#@ Integer (label = "Shrink nuclei (nr of erosions)", style = "spinner", min=0, max=100, value=2) erosions
#@ Float (label = "Minimum circularity", style = "spinner", min=0, max=1, value=0.33) minCircularity
#@ Boolean (label = "Manually delete/edit segmented nuclei?", value=false) editNuclei
#@ Integer (label = "Measure up to which percentile?", style = "scroll bar", min=0, max=100, value=25) PercentileValue

// Initialize
saveSettings;
initialize();
setBatchMode(true);
start = getTime();

processFile(input_file, output);

//Save Image and Results
//saveAs("Tif", output + File.separator + file + "_lifetime");
//selectWindow("Results");
//saveAs("results",output + File.separator + "Results_all_cells.csv");


end = getTime();
print("-------------------------------------------------------------------");
print("Finished processing "+roiManager("count")+" cells in "+d2s((end-start)/1000,1)+" seconds. ("+d2s(roiManager("count")/((end-start)/1000),1)+" cells per second)");

restoreSettings;



function processFile(input_file, output) {
	filename = File.getName(input_file);
	dir = File.getParent(input_file);
	Ext.setId(input_file);
	Ext.getSeriesCount(nr_series);
	//fileName = Ext.getCurrentFile(input_file);
	histogramTable = "Histograms of "+filename;
	Table.create(histogramTable);
//	Plot.create("Histograms", "gray value", "counts");
//	Plot.show();

	for(s=0 ; s<nr_series ; s++) {
		//if(nImages>0) run("Close All");
		Ext.setSeries(s);
		Ext.getSeriesName(seriesName);
		print("Processing file, series "+s+1+"/"+nr_series+" : "+input_file+"  ||  " + seriesName);
		run("Bio-Formats Importer", "open=["+input_file+"] autoscale color_mode=Grayscale view=Hyperstack stack_order=XYCZT series_" + s + 1);	//Wrong bookkeeping of BioFormats
		image = getTitle();
		getDimensions(width, height, channels, slices, frames);
		run("Grays");
		
		segmentNuclei(image, nucleiChannel);
		measureSignal(image, seriesName, signalChannel, PercentileValue, nResults);
	}
	Plot.create("Mean of lowest " + PercentileValue + " percentile", "x", "Mean of lowest " + PercentileValue + " percentile");
	Plot.add("Circle", Table.getColumn("Mean", "Results"));
	Plot.setStyle(0, "blue,#a0a0ff,1.0,Circle");
	Plot.show();


	Plot.create("Histograms","gray value","counts");
	selectWindow(histogramTable);
	columnNames = split(Table.headings,'\t');
	for(i=0; i<nResults; i++) {
		data = Table.getColumn(columnNames[i],histogramTable);
		//Generate line with random color
		color1 = toHex(random*255);
		color2 = toHex(random*255);
		color3 = toHex(random*255);
		Plot.setColor("#"+color1+color2+color3);
		Plot.add("line", data);
	}
	setBatchMode(false);
	Plot.show();

//	Plot.create("Percentile Normalized", "x", "Mean of lower " + PercentileValue + " percentile normalized to Mean of nucleus");
//	Plot.add("Circle", Table.getColumn("PercentileNormalized", "Results"));
//	Plot.setStyle(0, "red,#ffa0a0,1.0,Circle");
//	Plot.show();
//	Plot.create("Histograms", "gray value", "counts");
//	Plot.add("line",counts);
//	Plot.update();

	selectWindow("Results");
	Table.renameColumn("Mean", "Meanlowest_"+PercentileValue+"%");
	updateResults();
	saveAs("Results", output + File.separator + filename + "_Results.tsv");

	selectWindow(histogramTable);
	saveAs("Results", output + File.separator + filename + "_Histograms.tsv");
}


function initialize() {
	run("Bio-Formats Macro Extensions");
	run("Conversions...", "scale");
	setBackgroundColor(0, 0, 0);
	run("Colors...", "foreground=white background=black selection=#007777");
	setOption("BlackBackground", true);	//This is the important one
	roiManager("Reset");
	print("\\Clear");
	run("Clear Results");
	if(nImages>0) run("Close All");
	if(isOpen("Summary")) close("Summary");
	run("Set Measurements...", "area mean standard min integrated median limit display redirect=None decimal=3");
	if(!File.exists(output)) {
		create = getBoolean("The specified output folder "+output+" does not exist. Create?");
		if(create==true) File.makeDirectory(output);		//create the output folder if it doesn't exist
		else exit;
	}
}


function segmentNuclei(image, channel) {
	roiManager("reset");
	selectWindow(image);
	Stack.setChannel(channel);
	run("Duplicate...", "title=nuclei duplicate channels=" +channel);
	//run("Median 3D...", "x=" + median3Dradius + " y=" + median3Dradius + " z=" + median3Dradius);
	//if(medianFilterRadius > 0) run("Median...", "radius="+medianFilterRadius);
	if(medianFilterRadius > 0) run("Gaussian Blur...", "sigma="+medianFilterRadius);
	setAutoThreshold("Li dark");
	if(autoThresholdNuclei == false) {
		setBatchMode("show");
		waitForUser("Adjust threshold if necessary");
	}
	run("Convert to Mask");
	rename("Mask");
	run("Fill Holes");

	for(i=0; i<erosions; i++) run("Erode", "slice");
	run("Watershed");
	for(i=0; i<erosions; i++) run("Erode", "slice");
	
	run("Analyze Particles...", " size=" + lowerDiameterLimit + "-" + upperDiameterLimit + " circularity=" + minCircularity + "-1.00 exclude include add");
	close("Mask");
	selectWindow(image);
	roiManager("Show All");
	if(editNuclei==true) editROIs(image);
}


function measureSignal(image, seriesName, channel, percentile, m) {
	selectWindow(image);
	setBatchMode("show");
	//setBatchMode("hide");
	Stack.setChannel(channel);
	
	pixelTable = "Pixel data "+image;
	Table.create(pixelTable);
	row=0;

	for (i = 0; i < roiManager("count"); i++) {
		roiManager("select",i);
		//measure the mean of the whole nucleus
		List.setMeasurements();
		meanSignalWholeNucleus = List.getValue("Mean");

		//Put the histogram of the nucleus in a Table
		getHistogram(values, counts, 256);
		sum = sumArray(counts);
		counts = divideArraybyScalar(counts, sum);
		Table.setColumn(seriesName+" cell_"+i+1, counts, histogramTable);
		Table.update;
		
		//list pixel values
		getRawStatistics(nPixels, mean, min, max, std, histogram);
		for (j = 0; j < histogram.length; j++) {
			for (k = 0; k < histogram[j]; k++) {
				Table.set("cell_"+i+1, row, j, pixelTable);
				row++;
			}
		}
		Table.update;
		row=0;

		//Determine percentile threshold
		total = 0;
		t=0;
		while (total < nPixels*percentile/100) {
			total += histogram[t];
			t++;
		}
		setThreshold(0,t-1);
/*		
		List.setMeasurements();
		area = List.getValue("Area");
		mean = List.getValue("Mean");
		stddev = List.getValue("StdDev");
		minimum = List.getValue("Min");
		maximum = List.getValue("Max");
		intDen = List.getValue("IntDen");

		setResult("Mean", i, mean);
*/
		run("Measure");
		resetThreshold;
		setResult("Nucleus", m, i+1);
		//setResult("Var/Mean", m, pow(getResult("StdDev", m),2)/(getResult("Mean", m)) );
		//setResult("Min/Mean", m, getResult("Min", m)/getResult("Mean", m) );
		setResult("MeanWholeNucleus", m, meanSignalWholeNucleus);
		setResult("PercentileNormalized", m, getResult("Mean", m)/meanSignalWholeNucleus);
		m++;

	}
	saveAs("Tif", output + File.separator + image + "_segmented");
	close();
	
	Table.save(output + File.separator + image + "_pixelData.tsv");
	close(pixelTable);
	updateResults();
}



function editROIs(image) {
shift=1;
ctrl=2; 
rightButton=4;
alt=8;
leftButton=16;
insideROI = 32;

flags=-1;

selectWindow(image);
roiManager("Show All without labels");
setOption("DisablePopupMenu", true);
setBatchMode(true);
resetMinAndMax();
Stack.setDisplayMode("composite");
Stack.setChannel(signalChannel);
run("Green");
setBatchMode("show");
color_ROIs();
print("\\Clear");
print("Delete, combine and draw new ROIs. \n- Left clicking while holding CTRL deletes a ROI.\n- Select multiple ROIs with shift-left mouse button and right-click to merge them into one ROI. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n");
//showMessage("Delete, combine and draw new ROIs. \n- Left clicking while holding CTRL deletes a ROI.\n- Select multiple ROIs with shift-left mouse button and right-click to merge them into one ROI. \n- Draw new ROIs using the magic wand or the freehand tool and press 't' to add. \n- Press space bar when finished editing.\n\nThis information is also printed in the log window.");
print("\nStarting editing "+roiManager("count")+" ROIs...");


setTool("freehand");
roiManager("Show All without labels");
setOption("DisablePopupMenu", true);

nROIs = roiManager("Count");


while(!isKeyDown("space")) {				//exit by pressing space bar
	getCursorLoc(x, y, z, flags);
	if(flags==17 || flags==18)	{	// (De)select multiple ROIs with shift-leftclick; delete ROI with rightclick
		for(i=0;i<roiManager("Count");i++) {
			roiManager("Select",i);
			if(Roi.contains(x, y)==true) {
			selected = Roi.getProperty("selected");
				//click to select a single ROI
				if(flags==17 && selected==false) {		//select ROI
					//print("selecting ROI "+i);
					Roi.setStrokeColor("red");
					Roi.setProperty("selected",true);
				}
				else if(flags==17 && selected==true) {	//deselect ROI
					//print("deselecting ROI "+i);
					Roi.setStrokeColor("cyan");
					//Roi.setFillColor("1900ffff");
					Roi.setProperty("selected",false);
				}
				else if(flags==18) {	//delete ROI
					roiManager("Delete");
					for(j=0;j<roiManager("Count");j++) {	//deselect all ROIs and rename
						roiManager("Select",j);
						roiManager("Rename", "ROI "+j);
					}
				}
			}
		}
		roiManager("Deselect");
		run("Select None");
		updateDisplay();
	}


	if(flags==4) {	//right button: combine selected ROIs
		selected_ROI_array = newArray(roiManager("Count"));	//create array with indices of selected ROIs
		j=0;
		for(i=0;i<roiManager("Count");i++) {
			roiManager("select",i);
			selected = Roi.getProperty("selected");
			if(selected==true) {
				selected_ROI_array[j] = i;
				j++;
				//print(j);
			}
		}
		//check if more than one ROI is selected. If yes, combine the selected ROIs and update the list
		selected_ROI_array = Array.trim(selected_ROI_array,j);
		//print(selected_ROI_array.length + " ROIs selected");
		if(selected_ROI_array.length > 1) {
			//print("combining ROIs:");
			//Array.print(selected_ROI_array);
			roiManager("Select",selected_ROI_array);
			roiManager("Combine");
			roiManager("Update");
//			for(i=1;i<selected_ROI_array.length;i++) {	
			to_delete_array = Array.copy(selected_ROI_array);								//selecting and deleting redundant ROIs
			to_delete_array = Array.slice(selected_ROI_array,1,selected_ROI_array.length);	//create array without the first element
				roiManager("Deselect");
				//print("deleting ROIs:");
				//Array.print(to_delete_array);
				roiManager("select", to_delete_array);
				roiManager("Delete");
			roiManager("Select",selected_ROI_array[0]);
			//print("repairing ROI "+selected_ROI_array[0]);
			run("Enlarge...", "enlarge=1 pixel");			//remove wall between ROIs by enlarging and shrinking with 1 pixel
			run("Enlarge...", "enlarge=-1 pixel");
			roiManager("Update");
			
			setKeyDown("none");
			
			color_ROIs();
		}
	}


	if(nROIs!=roiManager("Count")) {	//change in the number of ROIs 
		run("Select None");
		color_ROIs();
		nROIs = roiManager("Count");
	}

	else wait(50);
}	//end of while loop

//Deselect and rename all ROIs once more
color_ROIs();
//print("Done editing. "+roiManager("count")+" ROIs remain.");
}


function color_ROIs() {
	run("Remove Overlay");

	for(j=0;j<roiManager("Count");j++) {	//fill all ROIs
		roiManager("Select",j);
		roiManager("Rename", "ROI "+j+1);
		Roi.setProperty("selected",false);
		//Roi.setFillColor("1900ffff");	//10% cyan fill
	}
	roiManager("Deselect");
	if(roiManager("count")>0) run("From ROI Manager");	//Add overlay containing the ROI fill
	roiManager("Select All");
//	roiManager("Set Color", "cyan");
	roiManager("Deselect");
	roiManager("Show All without labels");
	updateDisplay();
}

//Returns the sum of all elements of an arrays, neglecting NaNs
function sumArray(array) {
	sum=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
	}
	return sum;
}

//Divides all elements of an array by a scalar
function divideArraybyScalar(array, scalar) {
	divided_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		divided_array[a]=array[a]/scalar;
	}
	return divided_array;
}