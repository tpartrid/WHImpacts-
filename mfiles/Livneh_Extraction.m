%Livneh Extraction - based off script for NE Census Project 
% Extract Livneh Data to .mat file 
%   PART 1 - Extracts data from Livneh netcdf files (meteorological variables)
%   PART 2 - Extracts data from Livneh ascii files (sub daily fluxes)
%
%   **need to maually change script depending on variable  of interest**
clear; close all

netcdf = 0; %Run netcdf extraction
ascii = 1;  %Run ascii extraction 
sav = 1; %save output ? 

tic 

if netcdf == 1
    % Find Livneh data files
    direc='/adata/group/Livneh/FluxVars/';
    files=dir('/adata/group/Livneh/FluxVars/VIC_fluxes_Livneh_CONUSExt_v.1.2_2013.*.nc');

    % Read the lat, lon information for first year only
    % flipud lat because read in upside down relative to dataset
    ncid=netcdf.open(fullfile(direc,files(1).name),'NOWRITE');
    lon=netcdf.getVar(ncid,0);
    lat=flipud(netcdf.getVar(ncid,1));
    netcdf.close(ncid);

    %Create datevector based on desired start and end dates 
    %Livneh: daily data from 1915 01 01 - 2011 12 31 
    %Stored in monthly netcdf files 
    t1 = datenum(1915,1,1,0,0,0); %Make sure for loop below matches start date! 
    t2 = datenum(2015,12,31,0,0,0);
    t = cat(2,(t1:t2)',datevec(t1:t2));

    % Load preciptation and temp data into arrays. Preallocate size
    % RadIn_daily = zeros(444,922,length(t)); %Livneh is 444 x 922
    %Tx_daily = zeros(614,928,length(t));
    %Tn_daily = zeros(614,928,length(t));

    %Define CDL lat and long from metadata 
    % cdl_lat = linspace(51.6051354911,22.9401355905,82490)';
    % cdl_lon = linspace(-127.887227966,-65.3454731825,179978);
    % 
    % %Trim Livneh data to match cdl 
    % rows = lat < min(cdl_lat) | lat > max(cdl_lat);
    % cols = lon < min(cdl_lon) | lon > max(cdl_lon);
    % lat(rows) = [];
    % lon(cols) = [];

    % Load preciptation and temp data into arrays. Preallocate size
    %Pr_daily = zeros(614,928,length(t)); %Livneh is 444 x 922
    Qnet_daily = NaN(length(lat),length(lon),length(t));
    %Tn_daily = zeros(614,928,length(t));

        %Transform grid to column vectors of lon and lat   
        % [X,Y] = meshgrid(lon(1:10),lat(1:10)); %create Livneh grid
        % [Xq,Yq] = meshgrid(cdl_lon(301:330),cdl_lat(31:60));
        [X,Y] = meshgrid(lon,lat); %create Livneh grid
    %     [Xq,Yq] = meshgrid(cdl_lon,cdl_lat);
    %     Xq = single(Xq);
    %     Yq = single(Yq);
       % gc(:,1) = Xq(:); %lon
       % gc(:,2) = Yq(:); %lat
       % In_CNTY = false(length(gc),length(S)); %determine Logical array size 


    toc

    for k=1:1164 % 1930-2011, 15 years = 180, 97 yrs=1164
        k;
        %toc   
     % dir() include 4 dimesions:name,date,bytes,isdir
        %files(k).name
        ncid=netcdf.open(fullfile(direc,files(k).name),'NOWRITE');
    %    PrID = netcdf.inqVarID(ncid,'Prec');
        ID = netcdf.inqVarID(ncid,'Qnet');
    %    TnID = netcdf.inqVarID(ncid,'Tmin');
     %   pr_temp=netcdf.getVar(ncid,PrID); %temporary variable for each month
        tx_temp=netcdf.getVar(ncid,ID);
       % tn_temp=netcdf.getVar(ncid,TnID);

        %Find Index values to assign monthly values 
        Yr = str2double(files(k).name(end-8:end-5)); %retrieve year and month from file name
        M = str2double(files(k).name(end-4:end-3));
        fd = 1; %First day is always 1 
        ld = size(tx_temp,3); %Last day equal to length of temporary variable 
        idx_s = find(t(:,1) == datenum(Yr,M,fd)); %Compare datenums to index 
        idx_e = find(t(:,1) == datenum(Yr,M,ld));
      %  Pr_daily(:,:,idx_s:idx_e) = rot90(pr_temp);
       % Tn_daily(:,:,idx_s:idx_e) = rot90(tn_temp);
        tx_temp = rot90(tx_temp);
        size(tx_temp)
    %     tx_temp(rows,:,:) = [];
    %     tx_temp(:,cols,:) = [];
        tx_temp = single(tx_temp);
    %     size(tx_temp)    
    %     size(X)
    %     size(Y)
    %     size(Xq)
    %     size(Yq)

    %     Val = NaN(length(lat),length(lon),size(tx_temp,3));	
    % 	for d = 1:size(tx_temp,3)
    % 	    %Resample Livneh data to cdl resolution
    % % 	    IntVal = interp2(X,Y,squeeze(tx_temp(:,:,d)),Xq,Yq,'nearest');
    % %    	     CultVal = IntVal.*cdl;
    %     
    %     
    %    	 for j = 1:length(lat)
    % 		idxLat = cdl_lat <= lat(j) & cdl_lat > lat(j+1);
    % 		sum(idxLat)
    % 		for l = 1:length(lon)
    % 		    idxLon = cdl_lon >= lon(l) & cdl_lon < lon(l+1);
    % 	   	sum(idxLon) 
    % 		Valtemp = nanmean(CultVal(idxLat));
    % % 		Valtemp
    % 	  	  Val(j,l,d) = nanmean(Valtemp(idxLon));
    % 		end
    %   	 end
    %     end
        tx = tx_temp;
        Qnet_daily(:,:,idx_s:idx_e) = tx;
        netcdf.close(ncid);
    end

    % Replace missing values in raw data with NaNs
    %xp = find(Pr_daily(:,:,:)>1.0e19);
    xtx = find(Qnet_daily(:,:,:)>1.0e19);
    %xtn= find(Tn_daily(:,:,:)>1.0e19);
    %Pr_daily(xp) = NaN;
    %Tn_daily(xtn) = NaN;
    Qnet_daily(xtx) = NaN;



    % Save subsetted data
    if sav==1
       save /adata/tpartrid/Livneh_1915_2015/Qnet_latlon_Livneh.mat lat lon -v7.3
     %   save /adata/tpartrid/Livneh_1915_2015/Pr_Livneh.mat Pr_daily t -v7.3
     %   save /adata/tpartrid/Livneh_1915_2015/Tn_Livneh.mat Tn_daily t -v7.3
        save /adata/tpartrid/Livneh_1915_2015/Qnet_Livneh_.mat Qnet_daily t -v7.3
    end
end

%% Part 2 - Livneh ascii Extraction - Extracts sub daily Livneh data from ascii files 
% Extract Livneh Data to .mat file 

if ascii == 1
    load /adata/tpartrid/Livneh_1915_2015/QShortWave_dly.mat
    %Latitude and longitude arrasy 
    Y = (48.96875:-0.0625000:36.96875)';
    X = -104.96875:0.0625000:-79.96875;
    %Preset size of output arrays 
%     Qshrt_dly = NaN(length(Y),length(X),35429);
%     Qlong_dly = NaN(length(Y),length(X),35429); 
    format long 
    
    % Find Livneh data files
    % D=dir('/adata/group/Livneh/SubDaily/Run1/');
    % for d = 3:numel(D) %Loop through each folder in directory (ignore "." and "..")
    files=dir(('/adata/group/Livneh/SubDaily/Run1/fluxes.95.90.37.49/VIC_subdaily_fluxes_Livneh_CONUSExt_v.1.2_2013_43.03125_-90.78125'));
   % D(d).name9
    length(files)
    for i = 1:length(files)
        i;
        A = files(i).name; %Pull file name 
        lt = str2double(A(48:55)); %Pull latitude from file name
        ln = str2double(A(57:end)); %Pull longitude from file name
        m = dlmread(['/adata/group/Livneh/SubDaily/Run1/fluxes.95.90.37.49/',files(i).name]); %Read in file

        %Extract relevent variables from file
        QShrt = m(:,9).*10800; %3-hour shortwave in W/m^2 converted to J/m^2/3hours
   %     Qlng = m(:,11).*10800; %3-hour longwave in W/m^2 converted to J/m^2/3hours
        [a,~,c] = unique(m(:,1:3),'rows'); %Index days  IDing unique year / month / day combinations 
        ShrtDly = accumarray(c,QShrt,[],@nansum); %Total daily Q short in J/m^2/day
        ShrtDly = ShrtDly./1e6; %Convert to MJ
  %      LongDly = accumarray(c,Qlng,[],@nansum); %Total daily Q long in J/m^2/day
 %       LongDly = LongDly./1e6;


        %Put aggregated values in array
        y = find(Y==lt); x = find(X == ln);

        if ~isempty(y) && ~isempty(x)
            Qshrt_dly(y,x,:) = ShrtDly;
%            Qlong_dly(y,x,:) = LongDly; 
        end
    end
stop
    %Save updated files
    save /adata/tpartrid/Livneh_1915_2015/QShortWave_dly.mat Qshrt_dly -v7.3 %Save Daily incoming shortwave radiation, units: MJ/m^2/day
    %save /adata/tpartrid/Livneh_1915_2015/QLongWave_dly.mat Qlong_dly -v7.3 %Save Daily incoming longwave radiation, units: MJ/m^2/day
    toc
    stop
end
