%Aggregate Livneh data to county / ASD 
clear; close all;
tic  
%Livneh Files for interpolation 

% load /adata/tpartrid/CultivatedArea.mat
load /adata/tpartrid/Livneh_1915_2015/Scenarios/GDDKDD_grid_QMAP.mat
load /adata/tpartrid/Livneh_1915_2015/Scenarios/latlon_LivnehTrim.mat; %Livneh grid 
%  load /adata/tpartrid/Livneh_1915_2015/QShortWave_dly.mat
%  load /adata/tpartrid/Livneh_1915_2015/Orig_Livneh.mat;
% load AEI_CONUS.mat
toc
% %Load CDL Mask
%filename = '/adata/tpartrid/CDL/CDL_Bi_WGS84.tif'
 filename = '/adata/tpartrid/NASS/CDL_Bi_Mask.tif';

[CDLMask, CDL_R] = geotiffread(filename);
CDLMask = double(CDLMask);
CDLMask(CDLMask<1) = NaN;
CDL_lat = [CDL_R.LatitudeLimits(2) - 0.03125:-0.0625:CDL_R.LatitudeLimits(1) + 0.03125]';
CDL_lon = [CDL_R.LongitudeLimits(1) + 0.03125:0.0625:CDL_R.LongitudeLimits(2) - 0.03125];


%Define variables to interpolate
% rows = t(:,2) <= 1960;
val = GDD;
val2 = KDD;
clear GDD KDD; 

%Trim CDL Mask and Livneh data to match 

cols = CDL_lon < min(lon) | CDL_lon > max(lon);
rows = CDL_lat > max(lat) | CDL_lat < min(lat); %Livneh extends past CDL in latitude
CDL_lat(rows) = []; CDL_lon(cols) = []; 
CDLMask(rows,:) = []; CDLMask(:,cols) = [];% CDL_lat(end) = []; CDLMask(end,:) = [];

%val(rows,:,:) = [];
val = val.*CDLMask; %Mask out non cultivated areas 
val2 = val2.*CDLMask;
toc


%Load CDL Data mask for Linveh data 
%direc = '/adata/tpartrid/CDL/';
%files = dir('/adata/tpartrid/CDL/CultivatedArea.nc');
%ncid=netcdf.open(fullfile(direc,files(1).name),'NOWRITE');
%cdl = netcdf.getVar(ncid,2);
%netcdf.close(ncid);
%cdl = fliplr(rot90(cdl));
%toc
%Convert CDL to binary
%cdl(cdl<2) = NaN;
%cdl(cdl==2) = 1; 
%cdl = double(cdl);
%save /adata/tpartrid/CultivatedArea.mat cdl -v7.3

%Load shapefile
S = shaperead('/adata/tpartrid/NASS/Shapefiles/County/cb_2016_us_county_20m.shp');
%S = shaperead('/adata/tpartrid/CDL/Shapefiles/Export_Output.shp');
% S = shaperead('/adata/tpartrid/NASS/Shapefiles/AgDistrict/ASD_2012_20m.shp');


%Transform grid to column vectors of lon and lat  
%Try an part of data 
% [X,Y] = meshgrid(lon(1:10),lat(1:10)); %create Livneh grid
% [Xq,Yq] = meshgrid(cdl_lon(301:330),cdl_lat(31:60));
[X,Y] = meshgrid(lon,lat); %create Livneh grid
% [Xq,Yq] = meshgrid(cdl_lon,cdl_lat);
gc(:,1) = X(:); %lon
gc(:,2) = Y(:); %lat
In_CNTY = false(length(gc),length(S)); %determine Logical array size 
toc


% [pdsi_X,pdsi_Y] = meshgrid(pdsi_lon,pdsi_lat);
% for m = 1:size(pdsi,3)
% %     F = griddedInterpolant(pdsi_X, pdsi_Y, pdsi(:,:,m), 'cubic');
%     pdsi_q(:,:,m) = interp2(pdsi_X,pdsi_Y,pdsi(:,:,m),X_grid,Y_grid,'nearest');
% end


%Loop through counties, identifying grid cells in each 
for i = 1:length(S)
    i;
    X_CNTY = S(i).X;
    Y_CNTY = S(i).Y;
    In_CNTY(:,i) = inpolygon(gc(:,1),gc(:,2),X_CNTY',Y_CNTY');
end


%Find average Orig within each county 
%Reshape array to cell by time  
Valtemp = reshape(val,[length(gc),size(val,3)]);
Val2temp = reshape(val2,[length(gc),size(val2,3)]);
% PDSI = reshape(pdsi_q,[length(gc),size(pdsi,3)]);
% PDSI_ASD_mon = NaN(size(pdsi_q,3),length(S));
GDD_dly_QMAP = NaN(size(val,3),length(S));
KDD_dly_QMAP = NaN(size(val2,3),length(S));

% Delta_ASD_dly = NaN(size(val,3),length(S));


for c = 1:size(S,1) %Loop through counties
    %PDSI_ASD_mon(:,c) = nanmean(PDSI(In_CNTY(:,c),:));
    GDD_dly_QMAP(:,c) = nanmean(Valtemp(In_CNTY(:,c),:));
    KDD_dly_QMAP(:,c) = nanmean(Val2temp(In_CNTY(:,c),:));
 %    Delta_ASD_dly(:,c) = nanmean(Valtemp(In_CNTY(:,c),:));
end



save /adata/tpartrid/Livneh_1915_2015/County/Scenarios/GDDKDD_QMAP.mat GDD_dly_QMAP KDD_dly_QMAP t -v7.3
toc
% save /adata/tpartrid/AEI_ASD.mat Pr_ASD_CDL t -v7.3
STOP
%Check that it works 
%Set up map axes including projection, lat lon grid, and labels 
A = PDSI_ASD_mon;
%Load variable to plots into shape strucutre 
for k = 1:numel(S)
    S(k).Lon = S(k).X; %Convert X and Y to Lat and Lon for plotting 
    S(k).Lat = S(k).Y;
    %Load climate indices into structure 
    S(k).A = nanmean(A(:,k));
end

%Define Surface color for each 
colormap(parula) %(flipud(jet)) %colorcube 
surfaceColors = makesymbolspec('Polygon', {'A', ...
    [min([S.A]) max([S.A])], 'FaceColor', colormap}); %min([cut.mse]) max([cut.mse])

%Plot map 
scrsz = get(0,'ScreenSize');
hf = figure('Position',[scrsz(3)*.05 scrsz(4)*0.3 scrsz(3)*.95 scrsz(4)*.8]);
axesm('MapOrigojection','lambert','MapLonLimit',[-125 -67],'MapLatLimit',[25 49],'Grid','off','GColor','k','MLineLocation',2,'PLineLocation',2,...
        'MeridianLabel', 'off','ParallelLabel','off');
geoshow(S,'DisplayType', 'polygon', 'SymbolSpec', surfaceColors);
% tightmap
h = colorbar;
h.Label.String = 'Years';
h.Label.FontSize = 24;
h.FontSize = 22;
% set( h, 'YDir', 'reverse' );
caxis([min([S.A]) max([S.A])])
title('Length of Yield Record ','FontSize',26)
