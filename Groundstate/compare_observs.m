load("C:/data/groundstates_julia_grid_new_01_60_200000_10_96.0_working.mat")

% load("C:/data/groundstates_julia_grid_new_09_60_200000_60_98.0.mat")

E_julia  = Observ.EVec; 
mu_julia = Observ.mucVec;

load("C:/data/groundstates_T4_200000_60nK_large_055_60_200000_10_96.00.mat")
% load("C:/data/groundstates_T4_200000_60nK_large_117_60_200000_60_98.00.mat")
E_mat  = gather(Observ.EVec); 
mu_mat = gather(Observ.mucVec);

[E_julia, E_mat] = get_arrays(E_julia, E_mat);
[mu_julia, mu_mat] = get_arrays(mu_julia, mu_mat);

limits = 0;

figure; 

subplot(1, 2, 1)
hold on 
yyaxis left
plot(E_julia, "Color", "Red", "LineStyle", "-")
plot(E_mat, "Color", "Blue", "LineStyle", "-")
ylabel("Energy")
yyaxis right
plot(abs(E_julia - E_mat)./E_julia)
yline(0, "LineStyle", "--")
ylabel("Relative difference")
legend(["Julia", "Matlab", ""])
if limits == 1
    xlim([1, 5])
end

subplot(1, 2, 2)
hold on 
yyaxis left
plot(mu_julia, "Color", "Red", "LineStyle", "-")
plot(mu_mat, "Color", "Blue", "LineStyle", "-")
ylabel("Chemical potential")
yyaxis right
plot(abs(mu_julia - mu_mat)./mu_julia)
yline(0, "LineStyle", "--")
ylabel("Relative difference")
legend(["Julia", "Matlab", ""])
if limits == 1
    xlim([1, 5])
end