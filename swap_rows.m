function I = swap_rows(I,m,n)
% Function swaps row m of a matrix with row n

I([m n],:) = I([n m],:);

