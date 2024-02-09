scene = getArgument();
dir_path = "/home/mark/cellpose/HVN037a_BF_ALL_cellpose_masks/s" + scene + "_merged/";
img1_name = "HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene" + scene + "_t001_ORG_Merged.tif";
img1 = dir_path + img1_name;
open(img1);
open(img1);
run("Concatenate...", "open image1=img1_name image2=img1_name");

for (i = 2; i < 912; i++) {
        iWithZeros = IJ.pad(i, 3);
        img_name = "HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene" + scene + "_t" + iWithZeros + "_ORG_Merged.tif";
        img = dir_path + img_name;
        open(img);
        run("Concatenate...", "open image1=Untitled image2=img_name");
}
img_name = "HVN037a_NoStim_TFP_Venus_Cherry_IkBeKO_40XOil_RedBF-02_scene" + scene + "_t911_ORG_Merged.tif"
run("Concatenate...", "open image1=Untitled image2=img_name");
img_save = "/home/mark/Data/HVN037_BF_ALL_s" + scene + "_ORG_Concat.tif";
saveAs("Tiff", img_save);
close();

