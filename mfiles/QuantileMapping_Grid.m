%Quantile Mapping - Grid 
%This script alters Livneh climate observations through both quantile
%mapping and a delta change approach to remove the U.S. warming hole 
tic
clear; close all; 
path(path, '/ihome/tpartrid/Matlab/Toolboxes/'); % set path for functions 

precip = 0; %Set to 1 if adjusting precip 

% Load data - can likely only run one variable at a time due to size 
% load /adata/tpartrid/Livneh_1915_2015/Qnet_Livneh_.mat % *** IF MAPPING PR BE SURE TO CAHNGE qmap FUNCTION TO CONVERT NEGATIVE #'S TO "0" ***
load /adata/tpartrid/Livneh_1915_2015/QShortWave_dly.mat
% load /adata/tpartrid/Livneh_1915_2015/latlon_Livneh.mat   


toc
var = single(Qshrt_dly); %set variable 
clear Qshrt_dly; 
lat = Qshrt_lat; lon = Qshrt_lon; %Didn't generate Qshrt for entire ConUS, need to import Qshrt specific coords 

%Trim variable to domain (lat 37 -> 47 lon: -82 -> -102)
rows = lat < 36 | lat > 48;
cols = lon < -103 | lon > -81;
var(rows,:,:) = [];
var(:,cols,:) = [];
%Trim lat and lon
lat(rows) = [];
lon(cols) = [];

%Rename Orignal variable to save 
Qshrt_dly_Orig = var;
t = cat(2,datenum(t),t);
save /adata/tpartrid/Livneh_1915_2015/Scenarios/latlon_LivnehTrim.mat lat lon;
save /adata/tpartrid/Livneh_1915_2015/Scenarios/Qshrt_Orig.mat Qshrt_dly_Orig t -v7.3
% clear Qshrt_dly_Orig


% %Quantile Mapping 
Qshrt_dly_QMAP = NaN(size(var,1),size(var,2),size(var,3));
t_idx = t(:,2) <= 1957; %Index of rows priod to 1958 regime shift
for i = 1:size(var,1) %Loop thorugh rows
    i;
    for j = 1:size(var,2) %Loop through columns 
        j;
        REF = double(cat(2,t(t_idx,1),squeeze(var(i,j,t_idx)))); % Cat time and variable into reference distribution (distribution to be mapped to)
        ORIG = double(cat(2,t(t_idx ~= 1,1),squeeze(var(i,j,t_idx ~= 1)))); %Cat time and variable into original distribution (i.e. dist. to be mnapped)
        [Corrected, coeff] = QMAPP(REF,ORIG,precip); %call QMAPP function 
        Qshrt_dly_QMAP(i,j,t_idx~=1) = Corrected; %Rename mapped (i.e. corrected) data to variable of interest 
        Qshrt_dly_QMAP(i,j,t_idx) = var(i,j,t_idx);


        %2D Data - comment out outer for loop if doing 2D case 
%                 REF = double(cat(2,t(t_idx,1),squeeze(var(t_idx,j)))); % Cat time and variable into reference distribution (distribution to be mapped to)
%                 ORIG = double(cat(2,t(t_idx ~= 1,1),squeeze(var(t_idx ~= 1,j)))); %Cat time and variable into original distribution (i.e. dist. to be mnapped)
%                 [Corrected, coeff] = QMAPP(REF,ORIG); %call QMAPP function 
%                 Qshrt_dly_QMAP(t_idx~=1,j) = Corrected; %Rename mapped (i.e. corrected) data to variable of interest 
%                 Qshrt_dly_QMAP(t_idx,j) = var(t_idx,j);
    end
 end

Qshrt_dly_QMAP = single(Qshrt_dly_QMAP);
toc
save /adata/tpartrid/Livneh_1915_2015/Scenarios/Qshrt_QMAP.mat Qshrt_dly_QMAP t -v7.3
clear Qshrt_dly_QMAP
%Delta Change Method
%3D Data
var_PreShift = var(:,:,t_idx); % rename Qshrte 1957 shift
var_PostShift = var(:,:,t_idx ~=1); %Rename post 1957 shift 

%2D Data 
% var_QshrteShift = var(t_idx,:); % rename Qshrte 1957 shift
% var_PostShift = var(t_idx ~=1,:); %Rename post 1957 shift

sz = size(var_PostShift);
CorrectedVar = NaN(sz(1),sz(2),sz(3)); %Set variable to original values 
t_Pre = t(t_idx,:);
t_Post = t(t_idx~=1,:);

for i = 1:12 %Loop through months
    %3D Data 
    RefVarMon = var_PreShift(:,:,find(t_Pre(:,3) == i));
    ID_ORIGVarMon= find(t_Post(:,3) == i); 
    OrigVarMon = var_PostShift(:,:,ID_ORIGVarMon);
    
    %Delta adjustment, multiplicative for precip, additiive for temp and Qshrt 
    if precip == 1
        delta =nanmean(RefVarMon,3)./nanmean(OrigVarMon,3); 
        CorrectedVar(:,:,ID_ORIGVarMon) = OrigVarMon.*delta;
    else
        delta =nanmean(RefVarMon,3) - nanmean(OrigVarMon,3);
        CorrectedVar(:,:,ID_ORIGVarMon) = OrigVarMon+delta;
    end
    
%     %2D Data 
%     RefVarMon = var_QshrteShift(find(t_Qshrte(:,3) == i),:);
%     ID_ORIGVarMon= find(t_Post(:,3) == i); 
%     OrigVarMon = var_PostShift(ID_ORIGVarMon,:);
%     delta =nanmean(RefVarMon,1) - nanmean(OrigVarMon,1);
%     CorrectedVar(ID_ORIGVarMon,:) = OrigVarMon+delta;
    
    %Zero out any negative values (Precip and Qshort only)
    %CorrectedVar(CorrectedVar<0) = 0;
end
Qshrt_dly_Delta = single(cat(3, var_PreShift,CorrectedVar)); %Change to single to reduce file size 

% %Trim time to 1956-2015 
% t_idx = t(:,2) < 1956;
% t(t_idx,:) = [];
% Qshrt_dly_Orig(:,:,t_idx) = [];
% Qshrt_dly_QMAP(:,:,t_idx) = [];
% Qshrt_dly_Delta(:,:,t_idx) = [];
% toc

%Save Files 
save /adata/tpartrid/Livneh_1915_2015/Scenarios/Qshrt_Delta.mat Qshrt_dly_Delta t -v7.3
