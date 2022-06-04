%% Introduction to MATLAB programming
% Author: Ikuro Sato, Tokyo Tech, 4.26.2022

%% Exercise A
x1 = -180 : 1 : 180 ; % arithmetic progression
y1 = sin(x1 / 180 * pi) ; % sine function
figure(1)
plot(x1, y1)

%% Exercise B
x2 = 1 : 20 ;
y2 = exp(x2) ; % exponential
figure(2)
semilogy(x2, y2)

%% Exercise C
x3 = 1 : 100 ;
y3 = x3 .^ 2 ; % element-wise power
figure(3)
loglog(x3, y3)

%% Exercise D
x4 = [4, 6, 1, 3, 9, 8, 2, 7, 10, 5] ;
[val, idx] = max(x4) ;
disp([val, idx])
[val, idx] = min(x4) ;
disp([val, idx])
[val, idx] = sort(x4, 'ascend') ;
disp(val)
[val, idx] = sort(x4, 'descend') ;
disp(val)
mean(x4)
std(x4)
var(x4)
median(x4)

%% Exercise E
x5 = rand(2, 3) % random variables of 2x3 array sampled from Uniform Distribution
sum(x5, 1) % vertical sum
sum(x5, 2) % horizontal sum
x5(:) % vectorization
[2, 3] * x5 % row vector times 2x3-matrix
x5 * [2; 3; 4] % 2x3 matrix times column vector
x5 * [2, 3, 4]' % 2x3 matrix times transposed row vector
2 .* x5 % element-wise product
[2; 3] .* x5 % row-wise product
x5 .* [2, 3, 4] % column-wise product

%% Exercise F
x6 = magic(4) ; % 4x4 magic square
im2col(x6, [2, 2])
im2col(x6, [2, 3])
reshape(x6, 8, 2)

%% Exercise G
A = zeros(1, 10) ;
for k = 1 : 10
    if rem(k, 2) == 0
        A(k) = k^2 ;
    end
end
A

