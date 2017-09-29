# ParticleComplexity
A set of matlab scripts and functions to analyse microscopy images of aerosol particles. Basically a glorified wrapper for Matlab's Image Processing Toolbox

This script allows the user to select a .tif image (filename), then extracts all particles
within the image. For each particle, the area and perimeter are
calculated. Each particle is analyzed for distinct regions within the
particle. The number of interior regions, as well as the sum of the
interior region area and perimeter are calculated. The script has three
outputs, which appear in the directory where the script is run. 

1) filename_segmentation.pdf shows the results of the initial particle
segmentation
2) filename.txt is a tab-delimited text file with columns described in
the header info. Area and perimeter values are in units of pixels. The
first column gives the assigned number of each particle in the image for
matching with output images. 
3) folder named filename with segmented images of each particle to check
effectiveness of automated calculations

Dependencies: Matlab Image Processing Toolbox, export_fig
(http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
