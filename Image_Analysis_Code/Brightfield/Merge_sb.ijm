arg = split(getArgument(), "+");
scene = arg[0];
print(scene);
frame = arg[1];
print(frame);

// iWithZeros = IJ.pad(frame, 4);
// originalpath = "/home/brandon/Data/Mark/HVN032_BF_ALL_s" + scene + "_t0001-t3000/";
// orininalname = "HVN032_BF_ALL_s" + scene + "_t0001-t3000_t" + iWithZeros + "_ORG.tif";
// original = originalpath + orininalname;
// open(original);
// 
// maskpath = "/home/mark/cellpose/HVN032_BF_ALL_cellpose_masks/s" + scene + "/";
// maskname = "HVN032_BF_ALL_s" + scene + "_t0001-t3000_t" + iWithZeros + "_ORG_cp_masks.png";
// mask = maskpath + maskname;
// open(mask);
// run("Merge Channels...", "c4=" +  orininalname + " c1=" + maskname + " create");

// outputpath = "/home/mark/cellpose/HVN032_BF_ALL_cellpose_masks/s" + scene + "_merged/HVN032_BF_ALL_s" + scene + "_t0001-t3000_t" + iWithZeros + "_ORG_Merged.tif";
// saveAs("Tiff", outputpath);

// close();

