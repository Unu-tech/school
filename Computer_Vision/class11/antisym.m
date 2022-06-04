function matT = antisym(T)
%% antisymmetric matrix
% T is a 3-component vector. matT is the antisymmetric matrix of T.
matT = [
    0, -T(3), T(2) ;
    T(3), 0, -T(1) ;
    -T(2), T(1), 0] ;
end