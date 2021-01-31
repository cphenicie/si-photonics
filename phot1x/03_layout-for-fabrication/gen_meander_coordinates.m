% totalLength = 650 - (660.7 - 650);
totalLength = 450 - 10.7; % I'm not sure where this constant 10.7 comes from
xi = 14.6500;
yi = 129.59;
yf = -2.91000;
xmin = 40;
nlines = 10;
radCurve = 5;


extraToPad = xmin - xi;
narcs = 2 * (nlines - 1);
extraNeededPerArc = (2 * radCurve) - (2 * pi * radCurve / 4) ;

height = (yi - yf) / (nlines - 1);

width = (totalLength - 2 * extraToPad - height * (nlines - 1) + ...
         narcs * extraNeededPerArc) / nlines;
xfar = xmin + width;

coordStr = sprintf('%.5f\t%.5f\n%.5f\t%.5f\n', xi, yi, xfar, yi);

yline = yi - height;
for i = 1:nlines-2
    if mod(i, 2) == 0
        newStr = sprintf('%.5f\t%.5f\n%.5f\t%.5f\n', xmin, yline, xfar, yline);
        coordStr = [coordStr newStr];
    else
        newStr = sprintf('%.5f\t%.5f\n%.5f\t%.5f\n', xfar, yline, xmin, yline);
        coordStr = [coordStr newStr];
    end
    yline = yline - height;
end

finalStr = sprintf('%.5f\t%.5f\n%.5f\t%.5f\n', xfar, yline, xi, yline);
coordStr = [coordStr finalStr]

fid=fopen('coords.txt','w');
fprintf(fid, coordStr);
fclose(fid)
