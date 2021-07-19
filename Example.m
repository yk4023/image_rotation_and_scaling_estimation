%--------------------------------------------------------------------------
%  Examples of usage of the method proposed in:
%
%  [1] K. Yu, R.S. Yang, Z. Hui and A.J. Peng, "Joint estimation of image 
%  rotation angle and scaling factor,"
%
% The 'Lena' image used in this example can be downloaded from the USC-SIPI
% image database on:
% http://sipi.usc.edu/database/download.php?vol=misc&img=4.2.04
% For copyright information, please go to:
% http://sipi.usc.edu/database/copyright.php
%
%--------------------------------------------------------------------------
% This code is provided only for research purposes.
%--------------------------------------------------------------------------
% Clear all variables and close all figures
clear all; 
close all;
block_size = 256;
test_factor = 1.3;
test_angle = 40;

% get the gray image of Lena and produce the test image
image = imread('data/lena.tiff');
gray = rgb2gray(image);
rotation = imresize(gray, test_factor);
rot_sca = imrotate(rotation, test_angle);

% Obtaining cyclic correlation two-dimensional spectrum
Txx = CalculateTxx(rot_sca,block_size);
Txx = fftshift(Txx);

% High pass filtering for Txx
[h,w] = size(Txx);
Txx(h/2-3:h/2+3,w/2-3:w/2+3)=0;

% To (0,1)^2
Txx = fftshift(Txx);

% Get rotation angle and scaling factor
[angle,factor]=estimate_factor_and_angle(Txx);
fprintf('the picture ture angle is: %.2f,  pre angle is: %.2f\nture factor is: %.2f, factor is: %.2f\n', test_angle,angle,test_factor,factor);