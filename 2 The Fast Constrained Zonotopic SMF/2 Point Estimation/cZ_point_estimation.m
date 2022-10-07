function x = cZ_point_estimation(cZ)
%   Return a feasible point in a given constrained zonotope
%   (c) Yirui Cong, created: 12-Oct-2021, last modified: --


G_cZ = cZ.G;
c_cZ = cZ.c;
A_cZ = cZ.A;
b_cZ = cZ.b;
cwb_cZ = cZ.cwb;

xi = linprog([], [], [], A_cZ, b_cZ, -cwb_cZ', cwb_cZ');

x = G_cZ * xi + c_cZ;