% load("./compdata/initialize_matlab.h5"); 
% psi_mat = psi; 
% V_mat = V; 
% VDk_mat = VDk;
psi_mat     = h5read('./compdata/initialize_matlab.h5', '/psi');
V_mat       = h5read('./compdata/initialize_matlab.h5', '/V');
VDk_mat     = h5read('./compdata/initialize_matlab.h5', '/VDk');
dx_mat  = h5read('./compdata/initialize_matlab.h5', '/Transf/dx');
dy_mat  = h5read('./compdata/initialize_matlab.h5', '/Transf/dy');
dz_mat  = h5read('./compdata/initialize_matlab.h5', '/Transf/dz');

psi_julia       = h5read('./compdata/initialize_julia.h5', '/psi');
V_julia         = h5read('./compdata/initialize_julia.h5', '/V');
VDk_julia       = h5read('./compdata/initialize_julia.h5', '/VDk');
dx_julia    = h5read('./compdata/initialize_julia.h5', '/dx');
dy_julia    = h5read('./compdata/initialize_julia.h5', '/dy');
dz_julia    = h5read('./compdata/initialize_julia.h5', '/dz');

n_mat   = abs(psi_mat).^2;
nxz_mat = squeeze(sum(n_mat*dy_mat,2));
nyz_mat = squeeze(sum(n_mat*dx_mat,1));
nxy_mat = squeeze(sum(n_mat*dz_mat,3));

n_julia   = abs(psi_julia).^2;
nxz_julia = squeeze(sum(n_julia*dy_julia,2));
nyz_julia = squeeze(sum(n_julia*dx_julia,1));
nxy_julia = squeeze(sum(n_julia*dz_julia,3));

% figure; 

% subplot(2, 3, 1)
% plotxz = pcolor(x,z,nxz');
% set(plotxz, 'EdgeColor', 'none');
% xlabel('$x$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');
% daspect([1 1 1])
% title("Density integrated along $y$")
% colorbar;
% colormap(inferno)
   
% subplot(2, 3, 2)
% plotyz = pcolor(y,z,nyz');
% set(plotyz, 'EdgeColor', 'none');
% xlabel('$y$ [$\mu$m]'); ylabel('$z$ [$\mu$m]');
% daspect([1 1 1])
% title("Density integrated along $x$")
% colorbar;

% subplot(2, 3, 3)
% plotxy = pcolor(x,y,nxy');
% set(plotxy, 'EdgeColor', 'none');
% xlabel('$x$ [$\mu$m]'); ylabel('$y$ [$\mu$m]'); 
% title("Density integrated along $z$")
% % xlim([-7, 7]);
% % ylim([-7, 7]);
% daspect([1 1 1]);
% colorbar;