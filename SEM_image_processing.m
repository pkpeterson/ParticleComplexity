%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script processes SEM images to ascertain image complexity. This script
% allows the user to select a .tif image (filename), then extracts all particles
% within the image. For each particle, the area and perimeter are
% calculated. Each particle is analyzed for distinct particles within the
% particle. The number of interior particles, as well as the sum of the
% interior particle area and perimeter are calculated. The script has three
% outputs, which appear in the directory where the script is run. 
% 1) filename_segmentation.pdf shows the results of the initial particle
% segmentation
% 2) filename.txt is a tab-delimited text file with columns described in
% the header info. Area and perimeter values are in units of pixels. The
% first column gives the assigned number of each particle in the image for
% matching with output images. 
% 3) folder named filename with segmented images of each particle to check
% effectivness of automated calculations
% 
% Dependencies: Image Processing Toolbox, export_fig
% (http://www.mathworks.com/matlabcentral/fileexchange/23629-export-fig)
% PK Peterson 29 May 2015 (pkpeter@umich.edu)
% Modified 14 April 2016 pkp, now skips analysis of individual particles
% that are too homogenous to calculate interior structure. Interior values
% are returned as NaN

function SEM_image_processing(filename)
%clears command window and used variables
clc
clear M Inner_Area Inner_Perimeter num_interior_particles Outer_Area Outer_Perimeter
% opens file selection window, filters for .tif files, if you have other
% format images, this will need to change
% Uncomment if not using as a function
%filename=uigetfile('.tif');
I=imread(filename);
% Flatten 24 bit images
[x y z]=size(I);
if z>1
    I=rgb2gray(I);
end
% cuts off info bar at the bottom of each image, assumes height of 75 pixels
% if you know the actual height you should probably change this number
I=I(1:x-75,:);
% calculates thresholds for segmentation using Otsu's method (Otsu 1979,
% IEEE). I found using multiple level thresholding lead to better background
% separation than graythresh()
[level em]=multithresh(I,2)

% Performs segmentation of particles using threshold calculated above
seg_I=imquantize(I,level);
valuesMax = [level max(I(:))];
[quant8_I_max, index] = imquantize(I,level,valuesMax);
quant8_I_max(quant8_I_max>min(valuesMax))=max(valuesMax);
BW=im2bw(quant8_I_max);
%reduce noise for final particle count, eliminates particles smaller than 8
%pixels
cc = bwconncomp(BW);
stats = regionprops(cc, 'Area');
% this line sets threshold for minumum pixels in particle, so you can adapt
% it to your needs
idx = find([stats.Area] > 8);
BW = ismember(labelmatrix(cc), idx);
stats = regionprops(BW,'Area','Perimeter');
%get area and perimenter of all particles found
Outer_Area=[stats.Area];
Outer_Perimeter=[stats.Perimeter];


[L,num_particles]=bwlabel(BW);
if num_particles==0
    return
end
close all
%outputs image so you can verify particle segregation worked ok, figure is
%a pdf so it should be suitable for publications/presentations
tmp=strsplit(filename,'.');
filename=tmp{1};
output_image_name=strcat(filename,'_segmentation');
imshowpair(I,BW,'montage')
title('Image Segmentation')
export_fig(output_image_name, '-pdf')
%show number of particles in command window
num_particles


% Analysis of individual particles

% Makes folder for output images, if you have already run analysis on this
% image, this will overwrite previous analysis with current, it will throw
% a warning in the command window, but the script will still continue. If
% you feel inclined, you can modify this behavior. 
mkdir(filename);
cd(filename);

% Loops through all particles, calculates the number of interior particles,
% total area of interior particles, and total perimiter of interior
% particles
for i=1:num_particles
   close all
    mask=L;
    %Selects particle of interest
    mask(L~=i)=0;
    Outer_Area(i)
    Outer_Perimeter(i)
    mask(mask>0)=1;
    mask=uint8(mask);

    particle_image=I.*mask;
     
    orig=particle_image;
    % crops image
    [x y]=find(particle_image~=0);
    particle_image=particle_image(min(x):max(x),min(y):max(y));
    size(particle_image)
    orig=orig(min(x):max(x),min(y):max(y));
    % Enhances contrast for larger particles
    [a b]=size(particle_image);
    if a>8 && b>8
        particle_image = adapthisteq(particle_image);
    end

    
    
    % performs thresholding of individual particle
    [level em]=multithresh(particle_image,3)
    % Skips homogoneous particles as determined by the thresholds, the
    % number in the if statement be stricter (larger) or relaxed( smaller)
    % in integer values
    
    if (level(3)-level(2))<=45
             Inner_Area(i)=NaN;
             Inner_Perimeter(i)=0;
             num_interior_particles(i)=NaN;
             im_name=sprintf('Filtered%2.2d',i)
             imshow(orig);
             % comment line below if you don't feel an overwhelming need for
             % individual particle figures of homogeneous particles
             export_fig(im_name, '-png')
             continue
    end
    valuesMax = [level max(particle_image(:))];
    [quant8_particle_image, index] = imquantize(particle_image,level,valuesMax);






    % shows effectivness of thresholding. Interior particles are in white
    imshowpair(orig,quant8_particle_image,'montage')
    title('Particle Segmentation')
    quant8_particle_image(quant8_particle_image<max(valuesMax))=0;
    im_name=num2str(i);
    % comment line below if you don't feel an overwhelming need for
    % individual particle figures
    export_fig(im_name, '-png')
    %reduce noise for final particle count
    cc = bwconncomp(quant8_particle_image);
    stats = regionprops(cc, 'Area');
    idx = find([stats.Area] > 8);
    BW = ismember(labelmatrix(cc), idx);

    % calculates interior particles and properties of said particles
    [L_in,num_particles_in]=bwlabel(BW);

     num_particles_in
     num_interior_particles(i)=num_particles_in;
     stats_inner = regionprops(L_in,'Area','Perimeter');
     Inner_Area(i)=sum([stats_inner.Area]);
     Inner_Perimeter(i)=sum([stats_inner.Perimeter]);

end


% writes output file
M(:,1)=1:1:num_particles;
M(:,2)=Outer_Area;
M(:,3)=Outer_Perimeter;
M(:,4)=num_interior_particles;
M(:,5)=Inner_Area;
M(:,6)=Inner_Perimeter;

cd ../

outfile=strcat(filename,'.txt')
fid = fopen(outfile,'wt');
formatSpec='%s\t%s\t %s\t%s\t%s\t%s\n';
fprintf(fid, formatSpec,'Particle','Outer_Area','Outer_Perimeter','Num_Interior_Particles','Inner_Area','Inner_Perimeter');
formatSpec='%f\t%f\t %f\t%f\t%f\t%f\n';
[nrows,ncols] = size(M);
for row = 1:nrows
    fprintf(fid,formatSpec,M(row,:));
end
fclose(fid);
close all

end