%% Sample Code of Gaussian Filtering
%
% [Description] This program demostrates gaussian filtering.
%
% [Author] Ikuro Sato, Tokyo Tech, 3/14/2021

%% settings
r = -2 : 2 ; % range of filter (odd)
std_g = 2.0 ; % std for gaussian filter

%% gaussian kernel
[u, v] = meshgrid(r, r) ; % coordinates
tmp_f_GF = exp(-0.5 / std_g^2 * (u.^2 + v.^2)) ;
f_GF = tmp_f_GF ./ sum(tmp_f_GF(:)) ; % 2D gaussian filter

%% data loading
% Images may be found C:\Program Files\MATLAB\R2020b\toolbox\images\imdata
img = double(imresize(imread('llama.jpg'), 1/2)) / 255 ; % RGB (0-1 ranged)
[sy, sx, sc] = size(img) ; % image size
fs = numel(r) ; % filter size
colim_r = im2col(img(:,:,1), [fs, fs]) ; % colum-image (red)
colim_g = im2col(img(:,:,2), [fs, fs]) ; % colum-image (green)
colim_b = im2col(img(:,:,3), [fs, fs]) ; % colum-image (blue)

%% gaussian filtering: loop over spatial indices
img_GF_1 = zeros(sy-fs+1, sx-fs+1, sc) ; % initialization of gaussian-filtered image
tic
for i = 0 : sx - fs
    for j = 0 : sy - fs
        img_GF_1(j+1,i+1,1) = sum(sum(f_GF .* img(j+(1:fs), i+(1:fs), 1))) ; % gaussian-filtered image (red)
        img_GF_1(j+1,i+1,2) = sum(sum(f_GF .* img(j+(1:fs), i+(1:fs), 2))) ; % gaussian-filtered image (green)
        img_GF_1(j+1,i+1,3) = sum(sum(f_GF .* img(j+(1:fs), i+(1:fs), 3))) ; % gaussian-filtered image (blue)
    end
end
tpassed = toc ;
disp(['Loop version (s):  ', num2str(tpassed)])

%% gaussian filtering: no-loop
img_GF_2 = zeros(sy-fs+1, sx-fs+1, sc) ; % initialization of gaussian-filtered image
tic
img_GF_2(:,:,1) = reshape(f_GF(:)' * colim_r, sy-fs+1, sx-fs+1) ; % gaussian-filtered image (red)
img_GF_2(:,:,2) = reshape(f_GF(:)' * colim_g, sy-fs+1, sx-fs+1) ; % gaussian-filtered image (green)
img_GF_2(:,:,3) = reshape(f_GF(:)' * colim_b, sy-fs+1, sx-fs+1) ; % gaussian-filtered image (blue)
tpassed = toc ;
disp(['No-loop version (s):  ', num2str(tpassed)])

%% difference of two approach
dif = img_GF_2 - img_GF_1 ;
disp(['maximum difference:  ', num2str(max(abs(dif(:))))])

%% visualization
set(0, 'DefaultAxesFontSize', 15, 'DefaultAxesFontName', 'Yu Gothic')
figure(1)
imagesc(img)
title('Original Image')
figure(2)
imagesc(img_GF_2)
title('Gaussian-Filtered Image')
