%% Toy problem of 3D reconstruction
% [Description] This script demonstrates the essential matrix estimation, pose estimation,
% and 3D reconstruction from 2-view inputs. No outliers are included.
% [Author] Ikuro Sato, Tokyo Tech, 2022.

%% Settings
gt_O = 0.1 * randn(3, 1) ; % orientation (in rad) of the second camera wrt the first camera, groud truth
gt_T = 0.3 * randn(3, 1) ; % translation (in m), ground truth
NumPt = 32 ; % number of 3D points
FocLen = 100 ; % focal length (in pix)
cy = 0.0 ; % y-coordinate of the image center (in pix)
cx = 0.0 ; % x-coordinate of the image center (in pix)
Z_lower = 1.0 ; % Z-axis lower bound (in m) of datapoints
Z_scale = 4.0 ; % Z-axis scale (in m) of datapoints
XY_scale = 1.0 ; % X- or Y-axis scale (in m) of datapoints
obs_noise_scale = 0.1 * 0 ; % observation noise scale (in pix)

%% Data generation
IntMat = [FocLen, 0, cx ; 0, FocLen, cy ; 0, 0, 1] ; % internal matrix
Cam1 = IntMat * eye(3, 4) ; % projection of the 1st camera
Cam2 = IntMat * rotmat(gt_O)' * [eye(3), -gt_T] ; % projection of the 2nd camera
gt_X = [XY_scale * randn(2, NumPt) ; Z_lower + Z_scale * rand(1, NumPt) ; ones(1, NumPt)] ; % 3D datapoints
if 0 < sum(gt_X(3,:) < 0) + sum(Cam2(3,:) * gt_X < 0)
    keyboard
end
gt_x1 = [Cam1(1,:) * gt_X ./ (Cam1(3,:) * gt_X) ; Cam1(2,:) * gt_X ./ (Cam1(3,:) * gt_X) ; FocLen * ones(1, NumPt)] ; % ground-truth points that should be observed by camera 1
gt_x2 = [Cam2(1,:) * gt_X ./ (Cam2(3,:) * gt_X) ; Cam2(2,:) * gt_X ./ (Cam2(3,:) * gt_X) ; FocLen * ones(1, NumPt)] ; % ground-truth points that should be observed by camera 2
x1 = gt_x1 + [obs_noise_scale * randn(2, NumPt) ; zeros(1, NumPt)] ; % observation by camera 1
x2 = gt_x2 + [obs_noise_scale * randn(2, NumPt) ; zeros(1, NumPt)] ; % observation by camera 2
scale_x = std([x1(:) ; x2(:)]) ; % scale of x
nrm_x1 = x1 ./ scale_x ; % normalized observation
nrm_x2 = x2 ./ scale_x ; % normalized observation

%% Figure
set(0, 'defaultAxesFontSize', 16)
figure(1); clf; plot(gt_x1(1,:), gt_x1(2,:), '.', 'MarkerSize', 17); grid on; axis equal; set(gca, 'YDir', 'reverse'); title('Camea 1'); xlabel('x'); ylabel('y')
figure(2); clf; plot(gt_x2(1,:), gt_x2(2,:), '.', 'MarkerSize', 17); grid on; axis equal; set(gca, 'YDir', 'reverse'); title('Camea 2'); xlabel('x'); ylabel('y')

%% Essential matrix estimation
A = [nrm_x1 .* nrm_x2(1,:) ; nrm_x1 .* nrm_x2(2,:) ; nrm_x1 .* nrm_x2(3,:)] ; % coefficient matrix multiplied to essential matrix
[U_A, S_A, V_A] = svd(A, 'econ') ; % SVD on the coefficient matrix
EssentialMatrix = reshape(U_A(:, 9), 3, 3) ; % essential matrix: the last vector from U

%% Candidate pose estimation
[U_E, S_E, V_E] = svd(EssentialMatrix) ; % SVD on the essential matrix
est_T = U_E(:,3) ; % estimated T-hut: the 3rd vector of U up to scale (uts)
est_R{1} = det(U_E * rotmat([0, 0, pi/2]) * V_E') * U_E * rotmat([0, 0, -pi/2]) * V_E' ; % estimated rotational matrix
tmp = [U_E(:,2), U_E(:,1), U_E(:,3)] * rotmat([0,0,-pi/2]) * [V_E(:,2), V_E(:,1), V_E(:,3)]' ;
% This also works -> 
% tmp = U_E * rotmat([0,0,pi/2]) * V_E';
est_R{2} = det(tmp) * tmp ; % estimated rotational matrix

%% Triangulation
est_Cam2{1} = IntMat * est_R{1}' * [eye(3), -est_T] ; % candidate-1 projection matrix for camera 2
est_Cam2{2} = IntMat * est_R{1}' * [eye(3), est_T] ; % candidate-2 projection matrix for camera 2
est_Cam2{3} = IntMat * est_R{2}' * [eye(3), -est_T] ; % candidate-3 projection matrix for camera 2
est_Cam2{4} = IntMat * est_R{2}' * [eye(3), est_T] ; % candidate-4 projection matrix for camera 2
for hypo = 1 : 4 % loop over hypothesis (of signs)
    est_X{hypo} = zeros(4, NumPt) ; % reconstructed 3D points with candidate 1
    for k = 1 : NumPt % loop over data index
        B = [
            Cam1(1,:) - x1(1,k) * Cam1(3,:) ;
            Cam1(2,:) - x1(2,k) * Cam1(3,:) ;
            est_Cam2{hypo}(1,:) - x2(1,k) * est_Cam2{hypo}(3,:) ;
            est_Cam2{hypo}(2,:) - x2(2,k) * est_Cam2{hypo}(3,:) ] ; % coefficient matrix for (X, Y, Z, 1)'
        [U_B, S_B, V_B] = svd(B) ; % SVD on the coefficient matrix
        est_X{hypo}(:, k) = V_B(:,4) / V_B(4,4) ;
    end
end

%% Choosing correct hypothesis
for hypo = 1 : 4 % loop over hypothesis (of signs)
    cnt_neg_Z_from1 = sum(est_X{hypo}(3,:) < 0) ; % count the number of data in front of camera 1
    cnt_neg_Z_from2 = sum(est_Cam2{hypo}(3,:) * est_X{hypo} < 0) ; % count the number of data in front of camera 2
    repro_x2 = [
        est_Cam2{hypo}(1,:) * est_X{hypo} ./ (est_Cam2{hypo}(3,:) * est_X{hypo}) ;
        est_Cam2{hypo}(2,:) * est_X{hypo} ./ (est_Cam2{hypo}(3,:) * est_X{hypo})] ; % reprojection of estimated 3D
    RMSE = sqrt(mean(sum((repro_x2 - x2(1:2,:)).^2))) ; % root-mean-square-error of the reprojection
    fprintf('HYP-%d,  # negative Z from camera 1 = %d, # negative Z from camera 2 = %d, RMSE = %f\n', hypo, cnt_neg_Z_from1, cnt_neg_Z_from2, RMSE)
end
[est_R{1}' * rotmat(gt_O), est_R{2}' * rotmat(gt_O)]
[est_T ./ gt_T, -est_T ./ gt_T]

