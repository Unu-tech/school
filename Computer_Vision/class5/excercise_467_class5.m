%% Sample Code of Bilateral Filtering
% [description] This is a program to demonstrate bilateral filter with a color image.
% [reference] C. Tomasi and R. Manduchi, "Bilateral filtering for gray and color images", 
% ICCV 1998.
% https://users.soe.ucsc.edu/~manduchi/Papers/ICCV98.pdf
% [author] Ikuro Sato, Tokyo Tech
% [last update] 4/26/2022

%% settings
r = -10 : 1 : 10 ; % range of integral
sd_g = 2.5 ; % standard deviation for gaussian filter
sd_r = 10.0 ; % standard deviation for range filter (~0.2 for RGB, ~10 for Lab)
img_file = 'peppers.png' ;
flag_lab = true ; % use of Lab color space
set(0, 'DefaultAxesFontSize', 16)

%% function handles
rgb = @(x, fl) fl * lab2rgb(x) + ~fl * x ; % true; Lab->RGB, false: RGB->RGB

%% preprocessing
ImgOrig = double(imread(img_file)) / 255 ; % load original image
if flag_lab, ImgOrig = rgb2lab(ImgOrig); end
[sy, sx, sc] = size(ImgOrig) ; % image size
sf = numel(r) ; % size of filter
CI1 = im2col(ImgOrig(:,:,1), [sf, sf]) ; % column-image, R (L)
CI2 = im2col(ImgOrig(:,:,2), [sf, sf]) ; % column-image, G (a)
CI3 = im2col(ImgOrig(:,:,3), [sf, sf]) ; % column-image, B (b)
figure(1)
imagesc(rgb(ImgOrig, flag_lab)); title('Original image'); axis image

%% warmup: box filtering
tmp1 = sum(CI1, 1) / sf^2 ;
tmp2 = sum(CI2, 1) / sf^2 ;
tmp3 = sum(CI3, 1) / sf^2 ;
ImgBox = zeros(sy-sf+1, sx-sf+1, sc) ; % init. of the output image
ImgBox(:,:,1) = reshape(tmp1, sy-sf+1, sx-sf+1) ; % box-filtered image, R (L)
ImgBox(:,:,2) = reshape(tmp2, sy-sf+1, sx-sf+1) ; % box-filtered image, G (a)
ImgBox(:,:,3) = reshape(tmp3, sy-sf+1, sx-sf+1) ; % box-filtered image, B (b)
figure(2)
imagesc(rgb(ImgBox, flag_lab)); title('Box-filtered image'); axis image

%% Gaussian filtering
[x, y] = meshgrid(r, r) ; % coodinates
tmpGaussianFilter = exp(-0.5 / sd_g^2 * (x.^2 + y.^2)) ; % Gaussian filter
GaussianFilter = tmpGaussianFilter / sum(tmpGaussianFilter(:)) ;
tmp1 = GaussianFilter(:)' * CI1 ;
tmp2 = GaussianFilter(:)' * CI2 ;
tmp3 = GaussianFilter(:)' * CI3 ;
ImgGauss = zeros(sy-sf+1, sx-sf+1, sc) ; % init. of the output image
ImgGauss(:,:,1) = reshape(tmp1, sy-sf+1, sx-sf+1) ; % Gaussian-filtered image, R (L)
ImgGauss(:,:,2) = reshape(tmp2, sy-sf+1, sx-sf+1) ; % Gaussian-filtered image, G (a)
ImgGauss(:,:,3) = reshape(tmp3, sy-sf+1, sx-sf+1) ; % Gaussian-filtered image, B (b)
figure(3)
imagesc(rgb(ImgGauss, flag_lab)); title('Gaussian-filtered image'); axis image

%% bilateral filtering
cen = 0.5 * (1 + sf^2) ; % index of the center pixel
RangeFilter = exp(-0.5/sd_r^2 * ((CI1-CI1(cen,:)).^2 + (CI2-CI2(cen,:)).^2 + (CI3-CI3(cen,:)).^2)) ; % range filter
tmpBilateralFilter = GaussianFilter(:) .* RangeFilter  ;
BilateralFilter =  tmpBilateralFilter ./ sum(tmpBilateralFilter) ; % bilateral filter
tmp1 = sum(BilateralFilter .* CI1, 1) ;
tmp2 = sum(BilateralFilter .* CI2, 1) ;
tmp3 = sum(BilateralFilter .* CI3, 1) ;
ImgBilateral = zeros(sy-sf+1, sx-sf+1, sc) ; % init. of the output image
ImgBilateral(:,:,1) = reshape(tmp1, sy-sf+1, sx-sf+1) ; % bilateral-filtered image, R (L)
ImgBilateral(:,:,2) = reshape(tmp2, sy-sf+1, sx-sf+1) ; % bilateral-filtered image, G (a)
ImgBilateral(:,:,3) = reshape(tmp3, sy-sf+1, sx-sf+1) ; % bilateral-filtered image, B (b)
figure(4)
imagesc(rgb(ImgBilateral, flag_lab)); title('Bilateral-filtered image'); axis image
% Observe:
% surf(reshape(BilateralFilter(:,(196-10)*(sy-20) + 95-10), sf, sf))
% surf(reshape(BilateralFilter(:,(200-10)*(sy-20) + 95-10), sf, sf))
