function [error] = F_error_2d(npoints,pos,points)

BoxScale = 100;
if npoints == 0
    error = nan;
end
error = BoxScale*norm(pos-points)/npoints;
end

