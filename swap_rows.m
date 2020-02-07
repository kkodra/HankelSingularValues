function I = swap_rows(I,m,n)
% Function swaps row m of a matrix with row n
% A method used to construct a standard singularly perturbed system.
I([m n],:) = I([n m],:);

