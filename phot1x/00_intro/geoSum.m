function [s] = geoSum(r, n)
% Sums the first n terms of a geometric sequence where
% where the ratio between successive terms is r
if r == 1
    s=n;
else
    s = (1-r^n)/(1-r);
end
end