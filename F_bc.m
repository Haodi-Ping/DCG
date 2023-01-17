function [A] = F_bc(d_li,d_lj,d_lk,d_ij,d_ik,d_jk)
%% sub triangles
ljk = [0 1 1 1;1 0 d_lj^2 d_lk^2;1 d_lj^2 0 d_jk^2;1 d_lk^2 d_jk^2 0];
lki = [0 1 1 1;1 0 d_lk^2 d_li^2;1 d_lk^2 0 d_ik^2;1 d_li^2 d_ik^2 0];
lij = [0 1 1 1;1 0 d_li^2 d_lj^2;1 d_li^2 0 d_ij^2;1 d_lj^2 d_ij^2 0];
%% entire triangle
ijk = [0 1 1 1;1 0 d_ij^2 d_ik^2;1 d_ij^2 0 d_jk^2;1 d_ik^2 d_jk^2 0];
%% areas of triangle
s_ljk = sqrt(abs(-det(ljk)/16));
s_lki = sqrt(abs(-det(lki)/16));
s_lij = sqrt(abs(-det(lij)/16));
s_ijk = sqrt(abs(-det(ijk)/16));
%% coefficient a
a_li = s_ljk/s_ijk;
a_lj = s_lki/s_ijk;
a_lk = s_lij/s_ijk;
a_li = abs(a_li);
a_lj = abs(a_lj);
a_lk = abs(a_lk);
area_bias = [
    a_li+a_lj+a_lk-1;
    a_li+a_lj-a_lk-1;
    a_li-a_lj+a_lk-1;
    a_li-a_lj-a_lk-1;
    -a_li+a_lj+a_lk-1;
    -a_li+a_lj-a_lk-1;
    -a_li-a_lj+a_lk-1;
];
area_bias = abs(area_bias);
[~, ind_min] = min(area_bias(:));
sign_pattern = [
    1,1,1;
    1,1,-1;
    1,-1,1;
    1,-1,-1;
    -1,1,1;
    -1,1,-1;
    -1,-1,1;
];
A = sign_pattern(ind_min,:);
A = A.*[a_li a_lj a_lk];

end

