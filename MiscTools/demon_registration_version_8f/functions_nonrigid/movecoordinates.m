function coordOut=movecoordinates(coord,Fx,Fy,Fz,mode)
% This function movepixels, will (backwards) translate the pixels 
% of an 2D/3D image according to x, y (and z) translation images 
% (bilinear interpolated).
% The function is a wrapper around mex files movepixels_2d_double.c and
% movepixels_3d_double.c and movepixels_3d_single.c
%
% J = movecoordinates(coord,Tx,Ty,[],mode);
% 	or
% J = movecoordinates(coord,Tx,Ty,Tz,mode);
%
% Inputs;
%   Tx, Ty, Tz : The transformation images, describing the
%             (backwards) translation of every pixel in x,y and z direction.
%   mode: If 0: linear interpolation and outside pixels set to nearest pixel
%            1: linear interpolation and outside pixels set to zero
%            2: cubic interpolation and outsite pixels set to nearest pixel
%            3: cubic interpolation and outside pixels set to zero
%
%
% Outputs,
%   Iout : The transformed image
%
% Function is written by D.Kroon University of Twente (March 2009)

if(~exist('mode','var')), mode=0; end

if(size(coord,2)<3)
    roundCoordinate=round(coord);
    roundCoordinateIdx=sub2ind(size(Fx),roundCoordinate(:,2),roundCoordinate(:,1));
    deltaCoord=[Fx(roundCoordinateIdx),Fy(roundCoordinateIdx)];
    coordOut=coord+deltaCoord;
end



