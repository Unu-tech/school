%% Toy problem of 3D reconstruction
% [Description] This script demonstrates the factorization method of 3D reconstruction from 
% multiple affine camera images.  No outliers are included.
% The absolute position in 3D cannot be reconstructed.
% [Author] Ikuro Sato, Tokyo Tech, 2022.

%% Settings
scale_orientation = 0.1 ; % orientation scale (in rad)
scale_translation = 0.5 ; % translation scale (in m)
NumCam = 3 ; % number of cameras
NumPt = 32 ; % number of all 3D points
FocLen = 100 ; % focal length (in pix)
cy = 0.0 ; % y-coordinate of the image center (in pix)
cx = 0.0 ; % x-coordinate of the image center (in pix)
Z_lower = 2.0 ; % Z-axis lower bound (in m) of datapoints
obs_noise_scale = 0.001 ; % observation noise scale (in pix)

%% Data generation
IntMat = [FocLen, 0, cx ; 0, FocLen, cy] ; % internal matrix, same for all cameras
tmp = rand(3, NumPt) - [0.5, 0.5, 0.5]' ; % for unit cube generation
invfac = max(abs(tmp), [], 1) ;
unitcube = (invfac.^-1) .* tmp ; % unit cube (1m^3)
gt_X = [unitcube + [0, 0, 0.5 + Z_lower]' ; ones(1,NumPt)] ; % 3D points
Cam = cell(1, NumCam) ; % initialization of the camera matrices
Cam{1} = IntMat * eye(3, 4) ; % affine transformation by the 1st camera
gt_x = cell(1, NumCam) ; % initialization of the ground truth image positions
x = cell(1, NumCam) ; % initialization of the observation
gt_x{1} = Cam{1} * gt_X  ; % ground-truth points that should be observed by camera 1
for c = 2 : NumCam
    ori = scale_orientation * randn(3,1) ; % camera orientation wrt the 1st camera
    tra = scale_translation * randn(3,1) ; % camera translation wrt the 1st camera
    EucTrans = rotmat(ori)' * [eye(3), -tra] ; % Euclidean transformation from coord-1 to coord-c
    Cam{c} = IntMat * EucTrans ; % affine transformation by the camera c
    if 0 < sum(EucTrans(3,:) * gt_X < 0)
        keyboard
    end
    gt_x{c} = Cam{c} * gt_X ; % ground-truth points that should be observed by camera c
end
for c = 1 : NumCam
    x{c} = gt_x{c} + obs_noise_scale * randn(2, NumPt) ; % observation by camera c
end

%% Figure
for c = 1 : NumCam
    figure(10 + c); clf
    for i = 1 : NumPt
        plot(x{c}(1,i), x{c}(2,i), '.', 'MarkerSize', 15); hold on
    end
    grid on; axis equal; xlabel('x'); ylabel('y'); title(['Camera-', num2str(c)])
end

%% Preprocessing
std1 = std(x{1}(:)) ;
msnrm_x = cell(1, NumCam) ; % initialization of the mean-subtracted, normalized observation
for c = 1 : NumCam
    nrm_x = x{c} ./ std1 ; % normalized observation
    msnrm_x{c} = nrm_x - mean(nrm_x, 2) ; % mean-subtracted, normalized observation
end

%% Factorization method
D = zeros(2 * NumCam, NumPt) ; % data matrix
for c = 1 : NumCam
    D((1:2)+2*(c-1),:) = msnrm_x{c} ; 
end
[U, S, V] = svd(D, 'econ') ; % factorization
i = 1 : 2*NumCam ;
tmp_i = flipud(reshape(i, 2, NumCam)) ;
j = tmp_i(:)' ;
A = [
    U(:,1).*U(:,1), U(:,2).*U(:,1), U(:,3).*U(:,1), U(:,1).*U(:,2), U(:,2).*U(:,2), U(:,3).*U(:,2), U(:,1).*U(:,3), U(:,2).*U(:,3), U(:,3).*U(:,3) ;
    U(i,1).*U(j,1), U(i,2).*U(j,1), U(i,3).*U(j,1), U(i,1).*U(j,2), U(i,2).*U(j,2), U(i,3).*U(j,2), U(i,1).*U(j,3), U(i,2).*U(j,3), U(i,3).*U(j,3)] ;
B = [ones(2*NumCam, 1) ; zeros(2*NumCam, 1)] ;
QQ = pinv(A) * B ; % LS solution for QQ'
Q = chol(reshape(QQ, 3, 3))' ; % matrix Q
tmp_est_Cam = U(:,1:3) * Q ; % estimated camera matrix
third = cross(tmp_est_Cam(1,:), tmp_est_Cam(2,:)) ; % vertical component to the first two rows
K = [tmp_est_Cam(1:2,:) ; third] ;
est_Cam = tmp_est_Cam * K' ; % making the first camera identity
est_X = K * inv(Q) * S(1:3,1:3) * V(:,1:3)' ; % estimated 3d points

%% Confirmation
est_X ./ [gt_X(1:3,:)-mean(gt_X(1:3,:),2)]
figure(1)
    plot3(gt_X(1,:), gt_X(2,:), gt_X(3,:), '.')
    axis equal; grid on
    xlabel('X1')
    ylabel('X2')
    zlabel('X3')
    title('3D data')
figure(2)
    plot3(est_X(1,:), est_X(2,:), est_X(3,:), '.')
    axis equal; grid on
    xlabel('X1')
    ylabel('X2')
    zlabel('X3')
    title('3D reconstruction')

