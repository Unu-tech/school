%% Sample Code of Robust Estimation
% M-estimator of linear regression problem
% Ikuro Sato, Tokyo Tech, 5/24/2022

%% Ground truth of parameters
a0 = -0.5 ;
a1 = 1.0 ;
disp(['groud truth:  a0 = ', num2str(a0), ', a1 = ', num2str(a1)])

%% Data generation
N = 50 ; % # of inliers + # of outliers
N_out = 10 ; % # of outliers
fac_noise = 0.2 ; % noise factor for inliers
fac_out = 60.0 ; % noise factor for outliers
r = fac_noise * randn(N,1) ; % inlier noise
r(1:N_out) = fac_out * abs(r(1:N_out)) ; % outlier noise
x = sort(-20 + 40 * sort(rand(N,1))) ; % observation (x)
y = a0 + a1 * x + r ; % observation (y)
%----------------------------------------------------------------------
set(0, 'DefaultAxesFontSize', 17, 'DefaultAxesFontName', 'Yu Gothic UI')
figure(1); clf
plot(x, y, 'o')
grid on; axis equal
xlabel('x', 'FontSize', 17)
ylabel('y', 'FontSize', 17)
title('Robust line fitting', 'FontSize', 17)

%% Linear Least Squares
A = [ones(N,1), x] ; % data matrix
a_LS = inv(A'*A) * A' * y ; % linear least squares solution
disp(['Linear Least Squares:  a0 = ', num2str(a_LS(1)), ',  a1 = ', num2str(a_LS(2))])
line([min(x) max(x)], a_LS(1)+a_LS(2)*[min(x) max(x)], 'Color', 'g', 'LineWidth', 2)

%% IRLS (Bisquare a.k.a. Tukey)
c_bisquare = 4.685 ;
a_IRLS = a_LS ; % initialization of the solution by linear LS
disp('IRLS begins.')
for k = 1 : 15
    res = y - a_IRLS(1) - a_IRLS(2) * x  ; % residue
    res = res / (1.48 * median(abs(res - median(res)))) ; % rescaled residue
    w = (1 - (res/c_bisquare).^2).^2 ; % weight update
    w(abs(res)>c_bisquare) = 0 ;
    W = diag(w) ; % vector (w) -> diagonal matrix (W)
    a_IRLS = inv(A' * W * A + 1e-10 * eye(2)) * A' * W * y ; % weighted LS   
    disp(['k = ', num2str(k), ', a0 = ', num2str(a_IRLS(1)), ', a1 = ', num2str(a_IRLS(2))])
end
%----------------------------------------------------------------------
line([min(x) max(x)], a_IRLS(1)+a_IRLS(2)*[min(x) max(x)], 'Color', 'r', 'LineWidth', 2)
lh = legend('data', 'LS', 'IRLS', 'Location', 'SouthEast') ;
set(lh, 'FontSize', 13)
