%% Miniture image generation
% Ikuro Sato, Tokyo Tech, 2021.4.27

%% setting
filename = 'shibuya.JPG' ;
r = -7 : 1 : 7 ; % filter range
stdmax = 8.0 ; % maximum STD of the gaussian filter (pix)
stdmin = 0.1 ; % minimum STD of the gaussian filter (pix)
FacResize = 1/4 ; % factor of resizing
nrm = @(x) x./sum(x, 1) ; % vertical normalization
set(0, 'DefaultAxesFontSize', 16)

%% filtering
OrigImg = double(imresize(imread(filename), FacResize)) / 255.0 ; % load image
[sy, sx, sc] = size(OrigImg) ; % image size
[u, v] = meshgrid(r, r) ; % coodinates for the filter
[sv, su] = size(u) ; % filter size
step = 2 * (stdmax - stdmin) / (sy - sv) ;
sgm = abs(-stdmax+stdmin : step : stdmax-stdmin) + stdmin ; % standard deviation at y
GF = zeros(sv*su, sy-sv+1) ;
for k = 1 : sy - sv + 1 % loop over y
    GF(:,k) = nrm( exp(-0.5 / sgm(k)^2 * (u(:).^2 + v(:).^2)) ) ; % Gaussian filter with variable STD
end
W = repmat(GF, 1, sx-su+1) ; % weight matrix (in im2col form)
ci1 = sum(W .* im2col(OrigImg(:,:,1), size(u))) ;
ci2 = sum(W .* im2col(OrigImg(:,:,2), size(u))) ;
ci3 = sum(W .* im2col(OrigImg(:,:,3), size(u))) ;
MiniatureImg = zeros(sy-sv+1, sx-su+1, 3) ; % miniature image
MiniatureImg(:,:,1) = reshape(ci1, sy-sv+1, sx-su+1) ;
MiniatureImg(:,:,2) = reshape(ci2, sy-sv+1, sx-su+1) ;
MiniatureImg(:,:,3) = reshape(ci3, sy-sv+1, sx-su+1) ;

%% visualization
figure(1); clf
imagesc(OrigImg); axis image; title('Original')
figure(2); clf
imagesc(MiniatureImg); axis image; title('Fake miniature')
figure(3); clf
plot(1 : sy-sv+1, sgm); grid on; ylabel('STD (pix)'); xlabel('y'); title('STD(y)'); xlim([1, sy-sv+1])

