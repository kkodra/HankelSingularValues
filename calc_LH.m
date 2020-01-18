function [Li,Hi] = calc_LH(A1,A2,A3,A4,eps)
% Function computes matrices L and H needed for exact decoupling of SP system
% Newton method solution utilized for computation

A4inv = inv(A4);
L0 = A4inv*A3;

% Newton Method Solution of the L-equation

for i = 1:10
  Li = L0;
  D1i = A4+eps*Li*A2;
  D2i = -eps*(A1-A2*Li);
  Qi = A3+eps*Li*A2*Li;
  Li = lyap(D1i,D2i,-Qi);
  L0 = Li;
end
errNewton = eps*Li*A1-A4*Li-eps*Li*A2*Li+A3;

% H-equations via the direct method
Hi = lyap(D2i,D1i,-A2);

% Error test for H
HiError = Hi*A4-A2+eps*(Hi*Li*A2-A1*Hi+A2*Li*Hi);

