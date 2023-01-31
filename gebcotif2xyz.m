function gebcotif2xyz(fname1,fname2,flag1)
% ___________________________________
% This function gets the bathymetry data from
% a GEBCO geotiff file and the data are exported
% to two plain text files xyz. In one of them,
% the coordinates are geographical and the second,
% the coordinates are UTM (WGS84).
%
% GEBCO website:
% https://www.gebco.net/
% 
% inputs:
% fname1 = GEBCO geotiff filename.
% fname2 = xyz filename without extension.
% flag1 = flag to change the sign of the depth:
%         flag = 1 => depths are positive (using in Delft3d).
%
% Example:
%   gebcotif2xyz('baty.tiff','output',1);
% 
% This function requires Matlab version 2020 or more recent and
% Matlab mapping toolbox.
%
% Gabriel Ruiz
% IMTA-2022
% ver:1.00.2
%______________________________________________
narginchk(3,3);
% Knowing the Matlab version
v = version;
if str2num(v(18:21))>= 2020 && license('test','MAP_Toolbox') == 1
 
	% Getting the geospatial information from raster file
	[dArray,sRef] = readgeoraster(fname1,'OutputType','double');

    if isempty(sRef) == 0

	    % Getting the grid from the raster file
	    [latg,longg] = geographicGrid(sRef,"fullgrid");

	    % Transform the grid matrix to array with x y z
	    lon = longg(:);
	    lat = latg(:);
	    [x,y] = xyz_p_UTM(lat,lon);
	
	    if flag1 == 1
		    triadGeo = [lon lat -dArray(:)];
		    triadUtm = [x y -dArray(:)];
	    else
		    triadGeo = [lon lat dArray(:)];
		    triadUtm = [x y dArray(:)];
	    end
	    triadGeo = triadGeo';
	    triadUtm = triadUtm';

	    % Ploting
	    figure;
        pcolor(longg,latg,dArray);
	    shading flat;
	    colormap jet;
	    axis image;
	    title('Bathymetry from GEBCO (2D)','FontSize',14,'FontWeight','bold');
	    xlabel('Longitude, [°]','FontSize',12,'FontWeight','bold');
	    ylabel('Latitude, [°]','FontSize',12,'FontWeight','bold');
	    hc = colorbar;
	    set(get(hc,'label'),'string','depth [m]');
	    set(gca,'Box','on',...
			'XGrid','on',...
			'YGrid','on',...
			'FontSize',10,...
			'FontWeight','bold');
	    print('bathy2D.png','-dpng','-r600');

	    figure;
%     subplot(1,2,2)
	    surf(longg,latg,dArray);
	    shading flat;
	    colormap jet;
	    title('Bathymetry from GEBCO (3D)')
	    xlabel('Longitude, [°]');
	    ylabel('Latitude, [°]');
	    zlabel('Depth, [m]')
	    hc = colorbar;
	    set(get(hc,'label'),'string','depth [m]');
	    print('bathy3D.png','-dpng','-r600');

	    % Exporting the triad values
	    fid = fopen(horzcat(fname2,'GEO.xyz'),'w');
	    fid2 = fopen(horzcat(fname2,'UTM.xyz'),'w');
    
	    fprintf(fid,'%7.3f %7.3f %10.3f\r\n',triadGeo);
	    fprintf(fid2,'%10.3f %10.3f %10.3f\r\n',triadUtm);
	    fclose all;
    else
        error('This image is not a Geotiff image! Aborting...');
    end
else
	error('This Matlab version is not compatible to the script or Mapping toolbox is not loaded');
end
return

function [X,Y] = xyz_p_UTM(latitude,longitude)
% _______________________________________________________________
% With this little routine, you can get the transformation
% or the conversion of spherical coordinates to UTM 
% coordinates, using some WGS84.
% I'm using the Alberto Cotticia and Luciano Saruce's
% equations, these equations appear in  "Bolletino di
% Geodesia e Science Affini", num. 1.
%
% Inputs:
% latitude, longitude vectors.
%
% Outputs:
% X, Y.
%
% Author: Gabriel Ruiz Martinez, Civil Engineer.
%  Last Modification: May/05
% ______________________________________________________________
la=latitude;
lo=longitude;
sa = 6378137.000000; 
sb = 6356752.314245;
e2 = (((sa.^2)-(sb.^2)).^0.5)./sb;
e2cuadrada = e2.^2;
c = (sa.^2)/sb;
lat = la.*(pi/180);
lon = lo.*(pi/180);
Huso = fix((lo./6)+ 31);
S = ((Huso.*6)-183);
deltaS = lon-(S.*(pi/180));
a = cos(lat).*sin(deltaS);
epsilon = 0.5.*log((1+a)./(1-a));
nu = atan(tan(lat)./cos(deltaS))-lat;
v = (c./((1+(e2cuadrada.*(cos(lat)).^2))).^0.5).*0.9996;
ta = (e2cuadrada./2).*epsilon.^2.*(cos(lat)).^2;
a1 = sin(2.*lat);
a2 = a1.*(cos(lat)).^2;
j2 = lat+(a1./2);
j4 = ((3.*j2)+a2)./4;
j6 = ((5.*j4)+(a2.*(cos(lat)).^2))./3;
alfa = (3/4).*e2cuadrada;
beta = (5/3).*alfa.^2;
gama = (35/27).*alfa.^3;
Bm = 0.9996.*c.*(lat-alfa.*j2+beta.*j4-gama.*j6);
x = epsilon.*v.*(1+(ta./3))+500000;
y = nu.*v.*(1+ta)+Bm;
X = round(x.*1000)./1000;
Y = round(y.*1000)./1000;
return
