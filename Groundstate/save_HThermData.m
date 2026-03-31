load("Groundstate/HThermData.mat", "HThermFit");

coeffs_A = HThermFit.coeffFitType.A.params;
coeffs_B = HThermFit.coeffFitType.B.params;

filename = 'HThermFit_data.h5';

% Save coeffs_A
h5create(filename, '/coeffs_A', size(coeffs_A));
h5write(filename, '/coeffs_A', coeffs_A);

% Save coeffs_B
h5create(filename, '/coeffs_B', size(coeffs_B));
h5write(filename, '/coeffs_B', coeffs_B);



% X00 = coeffs_A(1); X10 = coeffs_A(2); X01 = coeffs_A(3); X20 = coeffs_A(4); X11 = coeffs_A(5); X02 = coeffs_A(6); 
% X1 = coeffs_B(1); X2 = coeffs_B(2); X3 = coeffs_B(3); X4 = coeffs_B(4); X5 = coeffs_B(5); X6 = coeffs_B(6); 
