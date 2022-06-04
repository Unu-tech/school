%% Sample Code of Robust Estimation
% M-estimator of non-linear regression problem
% Ikuro Sato, Tokyo Tech, 3/13/2021

%% plot settings
set(0, 'DefaultAxesFontSize', 17, 'DefaultAxesFontName', 'Yu Gothic UI')

%% Ground truth of parameters
a0 = -0.1 ;
a1 = -0.3 ;
a2 = 0.3 ;
a3 = 0.4 ;
disp(['groud truth:  a0 = ', num2str(a0), ', a1 = ', num2str(a1), ', a2 = ', num2str(a2), ', a3 = ', num2str(a3)])

%% Data generation
N = 60 ; % number of inliers + number of outliers
N_out = 4 ; % number of outliers
fac_noise = 0.03 ; % noise factor for inliers
fac_out = 60.0 ; % noise factor for outliers
r = fac_noise * randn(N,1) ; % inlier noises
r(1:N_out) = fac_out * abs(r(1:N_out)) ; % outlier noise
x = sort(-2 + 4 * sort(rand(N,1))) ; % observation (x)
y = a0 + a1*x + a2*x.^2 + a3*x.^3 + r ; % observation (y)
%----------------------------------------------------------------------
figure(2); clf
plot(x, y, 'o')
grid on; axis equal
xlabel('x', 'FontSize', 17)
ylabel('y', 'FontSize', 17)
title('Robust polynomial fitting', 'FontSize', 17)

%% Linear Least Squares
A = [ones(N,1), x, x.^2, x.^3] ; % data matrix
a_LS = inv(A'*A) * A' * y ; % linear least squares solution
disp(['Linear Least Squares:  a0 = ', num2str(a_LS(1)), ',  a1 = ', num2str(a_LS(2)), ', a2 = ', num2str(a_LS(3)), ', a3 = ', num2str(a_LS(4))])
%line([min(x) max(x)], a_LS(1)+a_LS(2)*[min(x) max(x)], 'Color', 'g', 'LineWidth', 2)
line(x, A*a_LS, 'Color', 'g', 'LineWidth', 2)

%% IRLS (Bisquare)
c_bisquare = 4.685 ;
a_IRLS = a_LS ; % ƒ‹[ƒv‚Ì‰Šú’l
disp('IRLS begins.')
for k = 1 : 15
    res = y - A * a_IRLS ; % residue
    res = res / (1.48 * median(abs(res - median(res)))) ; % rescaled residue
    w = (1 - (res/c_bisquare).^2).^2 ; % weight update
    w(abs(res)>c_bisquare) = 0 ;
    W = diag(w) ;
    a_IRLS = inv(A' * W * A) * A' * W * y ; % weighted LS    
    disp(['k = ', num2str(k), ', a0 = ', num2str(a_IRLS(1)), ', a1 = ', num2str(a_IRLS(2)), ', a2 = ', num2str(a_IRLS(3)), ', a3 = ', num2str(a_IRLS(4))])
end
%----------------------------------------------------------------------
line(x, A*a_IRLS, 'Color', 'r', 'LineWidth', 2)
lh = legend('data', 'LS', 'IRLS', 'Location', 'SouthEast') ;
set(lh, 'FontSize', 13)
