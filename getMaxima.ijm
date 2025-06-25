//looping over all images in a time series
// finding maxima to count mtDNA 
// jmc 01-04-24
//need to open which image file you're going to work on before running this code

t= getTitle();
count= getDirectory("folder for results");

//folder name
folder2 = count + File.separator + t  ;

//make a folder with this name with which to collect all subsequent datasets
File.makeDirectory(folder2);

run("Duplicate...", "duplicate channels=1");

//which prominence do we want this get maxima to run over, use the output of test below to set this. 
prominenC = 300;

for (i = 1; i <= nSlices; i++) {
    setSlice(i);
    // Prominence may need changing for each video!!
    run("Find Maxima...", "prominence=prominenC strict output=List");
 	saveAs("Results", folder2 +  File.separator + i + "-" + "MaximaCoordsP" + prominenC + ".csv");
	close("Results");
}

  
close();



// this code takes your image, runs Find maxima and prints the foci points over the top, saving the image so we can
// manually check the prominences in a given range, and choose the correct one. 
prom_start = 150;
prom_stop = 500;
increment = 50;
prominences = newArray(prom_start, prom_start+increment, prom_start+2*increment, prom_start+3*increment, prom_start+4*increment,prom_start+5*increment, prom_stop);

//small series of frames to test the prominences on. 
frameTests = newArray(1,50,100);

for (k = 0; k < lengthOf(prominences); k++) {

	
 for (g = 0; g < lengthOf(frameTests); g++) {
 	//define which frame and prominence value we're working with
 	
	prominence = prominences[k];
	f = frameTests[g];
	
	// grab the frame, and the channel (always channel 1 for SYBR)
	run("Duplicate...", "duplicate channels=1 frames="+f);
	current = getTitle();
	
	// Run find maxima- we're looking to get points output so we can compare to the real images, so we choose 'single points'
	// note that we are using the 'Strict' parameter
	run("Find Maxima...", "prominence=prominence strict output=[Single Points]");
	
    //must make points output 16 bit so we can merge it with microscope image
	run("16-bit");
    run("Merge Channels...", "c1=" + current + " c2=[" + current + " Maxima] create");
	
	//save
    saveAs("PNG", folder2  +  File.separator + "Composite" + f + "P-" + prominence + ".png");
	close();
 }
}





