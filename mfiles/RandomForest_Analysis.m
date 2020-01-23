%NASS Yield Analysis 
%   This script develops a random forest crop model and explores the imapct 
%   of the warming hole on corn yield.
%   Coutny level yield data aquired from USDA NASS
%
%
%The script is divided into multiple parts:
%       Part 1: Calculates planting and harvesting date based on historical
%       NASS data 
%       Part 2: Yield Climate Correlation Analysis (YieldClimCorr)
%           a. Finds a critical period by correlating Yield with climate
%           variables
%       Part 3: Data preparation:
%           a. calculates critical period and growing season averages and prepares data for
%           random forest 
%       Part 4: Random Forest to determine most important yield predictors
%           a. Input all climate indices into Random forest w/ county yield
%           as dependent variable
%      
%--------------------------------------------------------------------------
%clear; close all
path(path,'/ihome/tpartrid/Matlab/Toolboxes/')

%Switches
PlantingDate = 0; %Calculates average planting data for each county from historical data
YldClimCorr_Wkly = 0;
RF_DataPrep = 0; %Prepares daily climate data for input into random forest model
RandomForest = 1; %Runs random forest
    RunCounterFactual = 1; %Secondary Switch to run counterfactual scenarios through RF  

sav = 1; %Save output? 

%% County (state) planting date calculation 
%   Idnetify growing season start date (DOY) and growing season end (DOY)
%   GS_start = Last day of week in which >80% of state's corn is planted 
%   GS_End = Last Day of week in which >80% of state's corn is mature 
%   
%   Import NASS percent planted or percent mature data using import wizard.
%   Fields needed;
%       Year; Period; WeekEnding; State and State ANSI; Value (percent planted or mature) 
if PlantingDate == 1
    %Loop through states and years to return growing season start date 
    PlntDate = NaN((max(Year)-min(Year)+1),max(StateANSI)); %Year by StateANSI (col. # corresponds w/ state ANSI code) 
    for i = 1:max(StateANSI)
        st = find(StateANSI == i); %Rows associated with each state 
        y = Year(st); %Years each state has data 
        w = WeekEnding(st); %Weeks each state has data
        StateVal = Value(st);
        Abv50 = find(StateVal >= 50); %Find all values above 50% for all years 
        ystrt = min(y); yend = max(y);
        plntwk = w(Abv50);
        [C, ia, ic] = unique(y(Abv50)); %Return index correspondings with unique years (ia is first value >80% for a given year)
        idx = C - min(Year)+1;
        PlntDate(idx,i) = day(plntwk(ia),'dayofyear'); %Return day of year 
    end
    Year = unique(Year);

    %Save 
    save /adata/tpartrid/NASS/Corn_StMatureDate_19792017.mat PlntDate Year 
end

%% Weekly correlations with Yield 
if YldClimCorr_Wkly == 1
    load /adata/tpartrid/NASS/County/CGY_BUA_AllCNTY.mat
    load /adata/tpartrid/NASS/County/CGAH_CNTY_1930_2015.mat
    load /adata/tpartrid/NASS/Corn_StPlantDate_19792017.mat
    load /adata/tpartrid/NASS/Corn_StMatureDate_19812017.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Pr_Orig.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Tn_Orig.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Tx_Orig.mat
    S = shaperead('/adata/tpartrid/NASS/Shapefiles/County/cb_2016_us_county_20m.shp');
    
    %Rename climate variables and trim climate data to 1930-2011 
    rows = t(:,2)>=1930 & t(:,2) <= 2011; 
    Tn = Tn_dly_Orig(rows,:);
    Tx = Tx_dly_Orig(rows,:);
    Pr = Pr_dly_Orig(rows,:);
    yr = (1930:2011)';
    clear Tn_dly_Orig Tx_dly_Orig Pr_dly_Orig
    
    %Convert Daily climate vars to weekly
    PrWkAll = NaN(82,numel(S),53); TnuWkAll = NaN(82,numel(S),53); TxuWkAll = NaN(82,numel(S),53);
    TxxWkAll = NaN(82,numel(S),53); TnnWkAll = NaN(82,numel(S),53); gs_start = NaN(1,3220);
    gs_end = NaN(1,3220);
    
    wk = week(datetime(datevec(t(:,1)))); %Array Iding week for each DOY (i.e. 1-53)
    for i = 1930:2011 %Loops through each year 
        i;
        yridx = t(:,2)==i; %Index of rows in year 
        Prtemp = Pr(yridx,:); %Temporary yearly values 
        Tntemp = Tn(yridx,:);
        Txtemp = Tx(yridx,:);
        for w = 1:53 %Loop through 5 weeks bfore planting - 150 days after planting
            wkidx = wk(yridx)==w; % Return 365(6) x 1 week of year (i.e. 1-53) 
            PrWkAll(i-1930+1,:,w) = nansum(Prtemp(wkidx,:)); %Total weekly precip
            TnuWkAll(i-1930+1,:,w) = nanmean(Tntemp(wkidx,:)); %Average weekly Tn
            TxuWkAll(i-1930+1,:,w) = nanmean(Txtemp(wkidx,:)); %Average weekly Tx
            TxxWkAll(i-1930+1,:,w) = max(Txtemp(wkidx,:)); %Max weekly Tx 
            TnnWkAll(i-1930+1,:,w) = min(Tntemp(wkidx,:)); % Min weekly Tn 
        end
    end
    TaWkAll = (TxuWkAll + TnuWkAll)/2;
    
    %Convert planting date to gs start and end 
    for i = 2125%1:length(S)%Loops through each county
        state = str2double(S(i).STATEFP); %ID state county is in 
        if state <= 56 %All states have FIPS # < 56 
            gs_start(i) = ceil(nanmean(PlntDate(1:33,state)));% Percent planted data goes from 1979-2011
            gs_end(i) = ceil(nanmean(CornMatDate(1:31,state))); %Percent mature data goes from 1981-2011
        else 
            gs_start(i) = NaN; gs_end(i) = NaN; %S structure contains U.S. territories - set to NaN 
        end
    end
    %If there is no planting data, set growing season start and end to US average 
    gs_start(isnan(gs_start)) = 128; gs_end(isnan(gs_end)) = 227;
    %Convert to Week
    gs_start = ceil(gs_start./7); gs_end = ceil(gs_end./7);
    
    %Preset arrays for correlations analysis 
    Ta_wk1 = NaN(1,numel(S)); Ta_wk2 = NaN(1,numel(S)); Txx_wk1 = NaN(1,numel(S)); 
    Txx_wk2 = NaN(1,numel(S));Pr_wk1 = NaN(1,numel(S)); Pr_wk2 = NaN(1,numel(S));
    Prcorr_out = NaN(1,numel(S));Tacorr_out = NaN(1,numel(S)); Txxcorr_out = NaN(1,numel(S));

    %Loop through counties correlate detreded yield with detrended climate
    %variables for each week.
    for c = 2125%1:numel(S)
        c
        Y = NaN(82,1);
        Yield = CGY_AllCnty(1:82,c);
        md = ~isnan(Yield);
        Y(md) = detrend(Yield(md));
                
        %Pull weekly climate data for county specific growing
        %season from start of growing season to end of growing season 
        PrWk_Cnty = squeeze(PrWkAll(:,c,gs_start(c):gs_end(c)));  
        TaWk_Cnty = squeeze(TaWkAll(:,c,gs_start(c):gs_end(c)));
        TxxWk_Cnty = squeeze(TxxWkAll(:,c,gs_start(c):gs_end(c))); 

        %Make correlation table from 5 weeks before planting to 150 days after planting 
        for w1 = 1:size(TaWk_Cnty,2)
            for w2 = 1:size(TaWk_Cnty,2)
                ta = detrend(nanmean(TaWk_Cnty(:,w1:w2),2));
                txx = detrend(nanmean(TxxWk_Cnty(:,w1:w2),2));
                pr = detrend(nansum(PrWk_Cnty(:,w1:w2),2));

                Prcorr(w1,w2) = corr(Y,pr,'rows','complete');
                Tacorr(w1,w2) = corr(Y,ta,'rows','complete');
                Txxcorr(w1,w2) = corr(Y,txx,'rows','complete');
            end
        end
        clear pr tn tx w1 w2 
        %Define critical period as period of week(s) with highest correlation
        %Prcp 
        [mn1, wkidx] = max(abs(Prcorr));
        [~, Prwk2] = max(mn1);
        Pr_wk2(c) = Prwk2;
        Pr_wk1(c) = wkidx(Pr_wk2(c));
        Prcorr_out(c) = max(mn1);
        %Ta
        [mn1, wkidx] = max(abs(Tacorr)); 
        [~, Tawk2] = max(mn1);
        Ta_wk2(c) = Tawk2;
        Ta_wk1(c) = wkidx(Ta_wk2(c)); 
        Tacorr_out(c) = max(mn1);
        %Txx
        [mn1, wkidx] = max(abs(Txxcorr));
        [~, Txxwk2] = max(mn1);
        Txx_wk2(c) = Txxwk2;
        Txx_wk1(c) = wkidx(Txx_wk2(c));
        Txxcorr_out(c) = max(mn1);
    end 
    CritPeriod = cat(1, Ta_wk1, Ta_wk2, Txx_wk1, Txx_wk2, Pr_wk1, Pr_wk2);  
    CritPeriodVars = {'TaWk1';'TaWk2';'TxxWk1';'TxxWk2';'PrWk1';'PrWk2'}; 
stop
    save /adata/tpartrid/TBagger_ProcData/Input/County/YieldClimCorr_Corn_Week.mat CritPeriod CritPeriodVars Txxcorr_out Tacorr_out Prcorr_out

end

%% Climate Varaible calculation 
if RF_DataPrep == 1
    
    %Load Files
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Pr_QMAP.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Tn_Orig.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/Tx_Orig.mat
    load /adata/tpartrid/Livneh_1915_2015/County/Scenarios/GDDKDD_Orig.mat
    load /adata/tpartrid/NASS/Corn_StPlantDate_19792017.mat%Average planting date by state 
    load /adata/tpartrid/NASS/Corn_StMatureDate_19812017.mat%Average Maturation date by state 
    load /adata/tpartrid/NASS/County/CGY_BUA_AllCNTY.mat
    load /adata/tpartrid/NASS/County/CGAH_CNTY_1930_2015.mat
    load /adata/tpartrid/TBagger_ProcData/Input/County/YieldClimCorr_Corn_Week.mat
    load /adata/tpartrid/TBagger_ProcData/Input/County/AEI_CNTY.mat
    S = shaperead('/adata/tpartrid/CDL/Shapefiles/cb_2016_us_county_20m_WithCDLArea.shp'); %County Shapefile with CDL area added 

    %Rename variables and trim data to 1930-2011 
    Yield = CGY_AllCnty(1:82,:);
    rows = t(:,2)>=1930 & t(:,2) <= 2011; 
    Tn = Tn_dly_Orig(rows,:);
    Tx = Tx_dly_Orig(rows,:);
    Pr = Pr_dly_QMAP(rows,:);
    GDD = GDD_dly_Orig(rows,:);
    KDD = KDD_dly_Orig(rows,:);
    t = t(rows,:); %Time 
    clear Tn_dly_Orig  Tx_dly_Orig  Pr_dly_Orig GDD_dly_Orig  KDD_dly_Orig  rows%Remove old vars 
    
    % PART 1: Set domain, ID counties to include and calc. non-climatic variables.
    %       Non-Climatic variables: County latitude, County AEI Fraction
    CNTY_latmean = NaN(1,3220); InDomain = NaN(1,3220); %Preset array size 
    
    %Define Domain: 82 - 102W & 37 - 47N 
    boxlat = [37*ones(1,201) 37:.01:47 47*ones(1,201) 47:-.01:37]; %Latitude for North band
    boxlon = [-102:.1:-82 -82*ones(1,1001) -82:-.1:-102 -102*ones(1,1001)]; %Longitude 
    for i = 1:numel(S)
        lat = S(i).Y(~isnan(S(i).Y));
        lon = S(i).X(~isnan(S(i).X));
        [latmean,lonmean] = meanm(lat,lon); %Define centroid as geographic average 
        CNTY_latmean(i) = latmean;
        in = inpolygon(lonmean,latmean,boxlon,boxlat); %Is centroid in lat. band? 
        InDomain(i) = in; %logical array ID'ing districts in band        
    end
    clear in latmean lonmean lat lon boxlat boxlon %Clear old vars 
    
    %Determine counties to include in analysis.
    %   County must have > 30 years yield data, >10000 acres harvest, have climate data*, & be in domain
    %   *Some counties have yield data but no cultivated area (CDL) exclude them from analysis 
    IncludedCnty = InDomain==1 & nanmean(CGAH_CNTY)>10000 & sum(~isnan(Yield))>=30 & ~isnan(nanmean(Tx));

    %Calculate Irrigatedfraction by district and filter
    TotCropArea_CDL = [S.SUM]*900; %Sum - total cult. grids in each county (30mx30m) convert to m^2 
    TotCropArea_CDL = TotCropArea_CDL/10000; %Convert m^2 to ha 
    AEIFrac = ((nanmean(AEI_CNTY))./ TotCropArea_CDL) .* 100; 
%     for i = 1:numel(S)
%     S(i).Area = areaint(S(i).Y,S(i).X)*51007200000;%Area of District in ha earths surface = 51.01 e 9 ha  
%     if length(S(i).Area)~=1
%         i
%         S(i).Area = sum(S(i).Area);
%     end
%     S(i).AEIFrac = nanmean(AEI_CNTY(:,i))/S(i).Area;
%     end
%     AEIFrac = [S.AEIFrac];%(nanmean(AEI_CNTY)./ TotCropArea_CDL) .* 100; 
    AEIFrac(IncludedCnty ~= 1) = NaN; %Not perfect alingment between cdl and AEI - set counties with 
    
    %Calculate Day of Year (DOY) to for comparison with growing season 
    ds = datestr(datenum(t(:,2:4))); %Convert year, month, day to datestring
    dt = datetime(ds,'InputFormat','dd-MMM-yyyy'); %convert datestring to datetime 
    t(:,5) = day(dt,'dayofyear'); %Save DOY in t 
    clear ds dt 

    %PART II: Calculate growing season and critical period aggregated for climate variables
    %       Climate variables: Tx, Tn, Pr, GDD, KDD
    Tn_gs = NaN(size(Tn,1),length(S)); Tx_gs = NaN(size(Tx,1),length(S));Pr_gs = NaN(size(Pr,1),length(S)); 
    GDD_gs = NaN(size(GDD,1),length(S)); KDD_gs = NaN(size(KDD,1),length(S)); gs_start = NaN(1,3220);
    gs_end = NaN(1,3220); Tx_gsu = NaN(82,numel(S)); Tn_gsu = NaN(82,numel(S)); Pr_gsu = NaN(82,numel(S));
    GDD_tot = NaN(82,numel(S)); KDD_tot = NaN(82,numel(S)); TaWk = NaN(82,3220); TxxWk = NaN(82,3220); 
    TnuWk = NaN(82,3220); TnnWk = NaN(82,3220); PrWk = NaN(82,3220); TxuWk = NaN(82,3220);
    
    %Select growing season based on DOY when >80% of corn is planted
    %PlntDate array calculated in PlntDatePrep. 
    t_gs = NaN(size(t,1),length(S));
    for i = 1:length(S)%Loops through each county
        state = str2double(S(i).STATEFP); %ID state county is in 
        if state <= 56 %All states have FIPS # < 56 
            gs_start(i) = ceil(nanmean(PlntDate(1:33,state)));% col in PlntDate corresponds with state FIPS #
            gs_end(i) = ceil(nanmean(CornMatDate(1:31,state))); %Growing Season end  - Average last day of week with > 90% mature 
        else 
            gs_start(i) = NaN;
            gs_end(i) = NaN;
        end
    end
    gs_start(isnan(gs_start)) = 121; %If there is no planting data, set growing season start to May 1. 
    gs_end(isnan(gs_end)) = 277; %U.S. average 

    %Extract growing season days from full time series 
    for i = 1:numel(S)
        gs = t(:,5)>=gs_start(i) & t(:,5)<=gs_end(i); % ID DOY during growing season 
        if sum(gs) >0 
            Tn_gs(gs,i) = Tn(gs,i); %If county has gs data - pull corresponding climate data 
            Tx_gs(gs,i) = Tx(gs,i);
            Pr_gs(gs,i) = Pr(gs,i);
            t_gs(gs,i) = t(gs);
            GDD_gs(gs,i) = GDD(gs,i);
            KDD_gs(gs,i) = KDD(gs,i);
        end
    end; clear i state gs
    %Calculate Ta as average of Tx and Tn 
    Ta_gs = (Tx_gs+Tn_gs)/2;
    
    %Calculate growing season aggregates of Tx, Tn, Pr, GDD, and KDD 
    for j = 1:numel(S) %Loop through counties
        j;
        y = year(datetime(datevec(t_gs(:,j)))); %convert datenum to year 
        for i=1930:2011 %loop through years
            yridx = y == i;  %Average across year (already truncated to growing season months)
            Txtemp = Tx_gs(yridx,j); Tntemp = Tn_gs(yridx,j); Prtemp = Pr_gs(yridx,j);
            GDDtemp = GDD_gs(yridx,j); KDDtemp = KDD_gs(yridx,j);
            Tx_gsu(i-1930+1,j) = nanmean(Txtemp); %Average GS Tx
            Txx_gsx(i-1930+1,j) = max(Txtemp);
            Tn_gsu(i-1930+1,j) = nanmean(Tntemp); %Avergae GS Tn
            Pr_gsu(i-1930+1,j) = nansum(Prtemp); %Total GS Pr 
            GDD_tot(i-1930+1,j) = nansum(GDDtemp);
            KDD_tot(i-1930+1,j) = nansum(KDDtemp);
            clear KDDtemp GDDtemp Txtemp Tntemp Prtemp yridx
        end
    end
    %nansum returns 0 for all NaN - convert to NaN
    Pr_gsu(Pr_gsu==0) = NaN;  

    %Calculate critical period variables 
    cp = round(nanmean(CritPeriod(:,IncludedCnty),2))*7; %determine average critical period 
    for c = 1:numel(S) %loop through each county  
        for y = 1930:2011
            yridx = t(:,2) == y; %Index of yearly rows 
            Txtemp = Tx_gs(yridx,c); %Temporary variables 
            Tntemp = Tn_gs(yridx,c);
            Tatemp = Ta_gs(yridx,c);
            Prtemp = Pr_gs(yridx,c);

            TaWk(y-1930+1,c) = nanmean(Tatemp(gs_start(c)+cp(1):gs_start(c)+cp(2))); % Only ran correlation for variables included in RF - Ta, Txx and Pr 
            TxuWk(y-1930+1,c) = nanmean(Txtemp(gs_start(c)+cp(1):gs_start(c)+cp(2)));
            TxxWk(y-1930+1,c) = max(Txtemp(gs_start(c)+cp(3):gs_start(c)+cp(4)));
            TnuWk(y-1930+1,c) = nanmean(Tntemp(gs_start(c)+cp(1):gs_start(c)+cp(2))); 
            TnnWk(y-1930+1,c) = min(Tntemp(gs_start(c)+cp(1):gs_start(c)+cp(2)));
            PrWk(y-1930+1,c) = nansum(Prtemp(gs_start(c)+cp(5):gs_start(c)+cp(6)));
        end
    end
    
    %Set counties not inlcuded to NaN 
    CNTY_lat(:,IncludedCnty~=1) = NaN; AEIFrac(:,IncludedCnty~=1) = NaN; Tx_gsu(:,IncludedCnty~=1) = NaN; Tn_gsu(:,IncludedCnty~=1) = NaN; Pr_gsu(:,IncludedCnty~=1) = NaN;
    TxWk(:,IncludedCnty~=1) = NaN; TxxWk(:,IncludedCnty~=1) = NaN; PrWk(:,IncludedCnty~=1) = NaN; GDD_tot(:,IncludedCnty~=1) = NaN; KDD_tot(:,IncludedCnty~=1) = NaN;
    
  
    if sav == 1
        save /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_QMAPPR_HistTx.mat Txx_gsx gs_start gs_end IncludedCnty CNTY_latmean AEIFrac Tx_gsu Tn_gsu  Pr_gsu TaWk TxxWk PrWk t GDD_tot KDD_tot t_gs
        disp('File saved in /adata/tpartrid/TBagger_ProcData/Input/County/')
    end
end

%% Random Forest to Identify predictors 
if RandomForest ==1 
    path(path,'/ihome/tpartrid/Matlab/Toolboxes/')
    
    %Load Data 
    load /adata/tpartrid/NASS/County/CGY_BUA_AllCNTY.mat
    load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_Orig.mat    

    %Limit yield analysis to only counties w/ > 30000 acres of harvested
    Yield = CGY_AllCnty(1:82,:).*.0254.*2.47; %Convert yield from bu/acre  to MT/ha .
    Yield(:,IncludedCnty~=1) = NaN; 
    yr = (1930:2011)'; %Timespan
    clear CGY_AllCnty 

    %reshape arrays to 1 x m  - run random forest over all districts
    y = reshape(Yield,[],1);
    yr_temp = repmat(yr,1,3220); yr_temp(:,IncludedCnty~=1) = NaN; yr_in = reshape(yr_temp,[],1);
    CNTY_latmean_temp = repmat(CNTY_latmean,82,1); CNTY_latmean_temp(:,IncludedCnty~=1) = NaN; CNTY_latmean_in = reshape(CNTY_latmean_temp,[],1);
    AEI_temp = repmat(AEIFrac,82,1); AEI_temp(:,IncludedCnty~=1) = NaN; AEI_in = reshape(AEI_temp,[],1);
    %TxuWk_in = reshape(TxuWk(:,IncludedCnty),[],1);
   % TnuWk_in = reshape(TnuWk(:,IncludedCnty),[],1);
    TaWk_in = reshape(TaWk,[],1);
    PrWk_in = reshape(PrWk,[],1);
    TxxWk_in = reshape(TxxWk,[],1);
   % TnnWk_in = reshape(TnnWk(:,IncludedCnty),[],1);
    Tx_gsu_in = reshape(Tx_gsu,[],1);
    Tn_gsu_in = reshape(Tn_gsu,[],1);
    Pr_gsu_in = reshape(Pr_gsu,[],1); 
    GDD_in = reshape(GDD_tot,[],1);
    KDD_in = reshape(KDD_tot,[],1); 

    %Put all potential features in table     
    Xtable = table(yr_in, CNTY_latmean_in, AEI_in,TaWk_in,PrWk_in,TxxWk_in,Tx_gsu_in,Tn_gsu_in,Pr_gsu_in,GDD_in,KDD_in); %IrrFrac_in ,CDHD_in ,GDDCr_in,KDDCr_in
    Xtable.Properties.VariableNames = {'Year','CNTY_lat','AEI','TaWk','PrWk','TxxWk','TxGs','TnGs','PrGs','GDD','KDD',};%Set variable names in table 'IrrFrac' ,'GDDcr','KDDcr' ,'CDHD'
    IV = {'Year','CNTY_lat','AEI','TaWk','TxxWk','PrWk','PrGs','GDD','KDD'};
    DV = 'Yield';
    X_Hist = table2array(Xtable(:,IV)); %Pull IV's from table

%         %Bayseian Optimization
%         maxMinLS = 20;
%         minLS = optimizableVariable('minLS',[1,maxMinLS],'Type','integer');
%         numPTS = optimizableVariable('numPTS',[1,size(X,2)-1],'Type','integer');
%         hyperparametersRF = [minLS; numPTS];
%         %Return results 
%         results = bayesopt(@(params)oobErrRF(params,X,y),hyperparametersRF,...
%             'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);
%         bestOOBErr = results.MinObjective;
%         bestHyperparameters = results.XAtMinObjective;
    B = TreeBagger(700,X_Hist,y,'Method','regression','MinLeafSize',10,...
        'PredictorNames',IV,'ResponseName',DV,'OOBPredictorImportance','on',...
        'Surrogate','on','PredictorSelection','interaction-curvature');
   
    %TreeBagger Performance 
    idx = ~isnan(y); %ID missing values in yield 
    ModVal.Y = y; 
    ModVal.OOBPredict = NaN(length(y),1);
    ModVal.OOBPredict(idx) = oobPredict(B);
    P_imp = B.OOBPermutedPredictorDeltaError;
    yhat = oobPredict(B);
    SST = nansum((y(idx)-nanmean(y(idx))).^2);
    SSR = nansum((y(idx)-yhat).^2);
    r2 = max(0,1 - SSR/SST);
    r2

    %Partial Dependence plots for each variable: for RF model 
    set(gcf,'visible','off') %supress figures
    for xv = 1:length(IV)
        pt1 = linspace(quantile(X_Hist(:,xv),0.01),quantile(X_Hist(:,xv),0.99))';
        figure
        ax = plotPartialDependence(B,IV{xv},'QueryPoints',pt1);
        xData(xv,:) = ax.Children.XData;
        yData(xv,:) = ax.Children.YData;
        close; clear ax 
    end

    %ALE plots for each variable: for RF model 
    set(gcf,'visible','off') %supress figures
    for xv = 1:length(IV)
        xv
        figure
        ax = ALE(B,IV{xv},60);%,'QueryPoints',pt1);
        ALE_xData(xv,:) = ax.XData;
        ALE_yData(xv,:) = ax.YData;
        close all
    end
     
%        %ALE plots for each variable: for RF model 
%         figure
%         ax = ALE(B,{'TxxWk','AEI'},50);%,'QueryPoints',pt1);
%         ALE_xData(xv,:) = ax.XData;
%         ALE_yData(xv,:) = ax.YData;
%         close 

        % Plot 2 way PDPs
        %TxxWk and PrWk
        figure
        pt1 = linspace(min(TaWk_in),max(TaWk_in),50)';
        pt2 = linspace(min(PrWk_in),max(PrWk_in),50)';
        ax = plotPartialDependence(B,{'TaWk','PrWk'},'QueryPoints',[pt1 pt2]);
        TaWkPrWk_xData(:,:) = ax.Children.XData;
        TaWkPrWk_yData(:,:) = ax.Children.YData;
        TaWkPrWk_zData(:,:) = ax.Children.ZData;
        close 
%         
%         TxGs and PrGs
        figure 
        pt1 = linspace(min(GDD_in),max(GDD_in),50)';
        pt2 = linspace(min(Pr_gsu_in),max(Pr_gsu_in),50)';
        ax = plotPartialDependence(B,{'GDD','PrGs'},'QueryPoints',[pt1 pt2]);
        PrGs_GDD_xData(:,:) = ax.Children.XData;
        PrGs_GDD_yData(:,:) = ax.Children.YData;
        PrGs_GDD_zData(:,:) = ax.Children.ZData;
%         close 

        % GDD and AEI
        figure 
        pt1 = linspace(min(GDD_in),max(GDD_in),50)';
        pt2 = linspace(min(AEI_in),max(AEI_in),50)';
        ax = plotPartialDependence(B,{'GDD','AEI'},'QueryPoints',[pt1 pt2]);
        AEI_GDD_xData(:,:) = ax.Children.XData;
        AEI_GDD_yData(:,:) = ax.Children.YData;
        AEI_GDD_zData(:,:) = ax.Children.ZData;
        
        % PrGS and AEI
        figure 
        pt1 = linspace(min(Pr_gsu_in),max(Pr_gsu_in),50)';
        pt2 = linspace(min(AEI_in),max(AEI_in),50)';
        ax = plotPartialDependence(B,{'PrGs','AEI'},'QueryPoints',[pt1 pt2]);
        AEI_PrGs_xData(:,:) = ax.Children.XData;
        AEI_PrGs_yData(:,:) = ax.Children.YData;
        AEI_PrGs_zData(:,:) = ax.Children.ZData;

    if RunCounterFactual == 1          
        %DELTA Scenario
        load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_Delta.mat  
        %reshape arrays to mx1 - keep original non-climatic variabels static 
        TaWk_in = reshape(TaWk,[],1);
        PrWk_in = reshape(PrWk,[],1);
        TxxWk_in = reshape(TxxWk,[],1);
        Tx_gsu_in = reshape(Tx_gsu,[],1);
        Tn_gsu_in = reshape(Tn_gsu,[],1);
        Pr_gsu_in = reshape(Pr_gsu,[],1); 
        GDD_in = reshape(GDD_tot,[],1);
        KDD_in = reshape(KDD_tot,[],1); 
        %Put all potential features in table     
        Xtable = table(yr_in, CNTY_latmean_in, AEI_in,TaWk_in,PrWk_in,TxxWk_in,Tx_gsu_in,Tn_gsu_in,Pr_gsu_in,GDD_in,KDD_in); %IrrFrac_in ,CDHD_in ,GDDCr_in,KDDCr_in
        Xtable.Properties.VariableNames = {'Year','CNTY_lat','AEI','TaWk','PrWk','TxxWk','TxGs','TnGs','PrGs','GDD','KDD'};%Set variable names in table 'IrrFrac' ,'GDDcr','KDDcr' ,'CDHD'
        X_Delta = table2array(Xtable(:,IV)); %Pull IV's from table
        ModVal.AltY_Delta = predict(B,X_Delta);

        % SCENARIO
        load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_QMAP.mat
        %reshape arrays to 1 x m  - run random forest over all districts 
        TaWk_in = reshape(TaWk,[],1);
        PrWk_in = reshape(PrWk,[],1);
        TxxWk_in = reshape(TxxWk,[],1);
        Tx_gsu_in = reshape(Tx_gsu,[],1);
        Tn_gsu_in = reshape(Tn_gsu,[],1);
        Pr_gsu_in = reshape(Pr_gsu,[],1); 
        GDD_in = reshape(GDD_tot,[],1);
        KDD_in = reshape(KDD_tot,[],1); 
        %Put all potential features in table     
        Xtable = table(yr_in, CNTY_latmean_in, AEI_in,TaWk_in,PrWk_in,TxxWk_in,Tx_gsu_in,Tn_gsu_in,Pr_gsu_in,GDD_in,KDD_in); %IrrFrac_in ,CDHD_in ,GDDCr_in,KDDCr_in
        Xtable.Properties.VariableNames = {'Year','CNTY_lat','AEI','TaWk','PrWk','TxxWk','TxGs','TnGs','PrGs','GDD','KDD'};%Set variable names in table 'IrrFrac' ,'GDDcr','KDDcr' ,'CDHD'
        X_QMAP = table2array(Xtable(:,IV)); %Pull IV's from table
        ModVal.AltY_QMAP = predict(B,X_QMAP);
        
        %QMAP Tx Historical Pr 
        load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_QMAPTx_HistPr.mat
        %reshape arrays to mx1 - keep original non-climatic variabels static 
        TaWk_in = reshape(TaWk,[],1);
        PrWk_in = reshape(PrWk,[],1);
        TxxWk_in = reshape(TxxWk,[],1);
        Tx_gsu_in = reshape(Tx_gsu,[],1);
        Tn_gsu_in = reshape(Tn_gsu,[],1);
        Pr_gsu_in = reshape(Pr_gsu,[],1); 
        GDD_in = reshape(GDD_tot,[],1);
        KDD_in = reshape(KDD_tot,[],1); 
        %Put all potential features in table     
        Xtable = table(yr_in, CNTY_latmean_in, AEI_in,TaWk_in,PrWk_in,TxxWk_in,Tx_gsu_in,Tn_gsu_in,Pr_gsu_in,GDD_in,KDD_in); %IrrFrac_in ,CDHD_in ,GDDCr_in,KDDCr_in
        Xtable.Properties.VariableNames = {'Year','CNTY_lat','AEI','TaWk','PrWk','TxxWk','TxGs','TnGs','PrGs','GDD','KDD'};%Set variable names in table 'IrrFrac' ,'GDDcr','KDDcr' ,'CDHD'
        X_QMAPTx = table2array(Xtable(:,IV)); %Pull IV's from table
        ModVal.AltY_QMAPTx = predict(B,X_QMAPTx);
        
        %QMAP Pr Historical Tx 
        load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures_QMAPPR_HistTx.mat         
        %reshape arrays to mx1 - keep original non-climatic variabels static 
        TaWk_in = reshape(TaWk,[],1);
        PrWk_in = reshape(PrWk,[],1);
        TxxWk_in = reshape(TxxWk,[],1);
        Tx_gsu_in = reshape(Tx_gsu,[],1);
        Tn_gsu_in = reshape(Tn_gsu,[],1);
        Pr_gsu_in = reshape(Pr_gsu,[],1); 
        GDD_in = reshape(GDD_tot,[],1);
        KDD_in = reshape(KDD_tot,[],1); 
        %Put all potential features in table     
        Xtable = table(yr_in, CNTY_latmean_in, AEI_in,TaWk_in,PrWk_in,TxxWk_in,Tx_gsu_in,Tn_gsu_in,Pr_gsu_in,GDD_in,KDD_in); %IrrFrac_in ,CDHD_in ,GDDCr_in,KDDCr_in
        Xtable.Properties.VariableNames = {'Year','CNTY_lat','AEI','TaWk','PrWk','TxxWk','TxGs','TnGs','PrGs','GDD','KDD'};%Set variable names in table 'IrrFrac' ,'GDDcr','KDDcr' ,'CDHD'
        X_QMAPPr = table2array(Xtable(:,IV)); %Pull IV's from table
        ModVal.AltY_QMAPPr = predict(B,X_QMAPPr);
        
        %Hold year constant to 1930
        X_Hist_noyr = X_Hist;
        X_Hist_noyr(:,1) = 1930;
        ModVal.AltY_Notech = predict(B,X_Hist_noyr);
    end

    if sav == 1 
        save /adata/tpartrid/TBagger_ProcData/Output/TBOut_EntireDomain_trended.mat r2 P_imp X_Hist X_Delta X_QMAP xData...
            yData IV ModVal gs_start gs_end TaWkPrWk_xData TaWkPrWk_yData TaWkPrWk_zData  PrGs_GDD_xData PrGs_GDD_yData PrGs_GDD_zData... 
            AEI_GDD_xData AEI_GDD_yData AEI_GDD_zData AEI_PrGs_xData AEI_PrGs_yData AEI_PrGs_zData ALE_xData ALE_yData -v7.3 % TaWkPrWk_xData TaWkPrWk_yData TaWkPrWk_zData  PrGs_GDD_xData PrGs_GDD_yData PrGs_GDD_zData
    end
end





 
  

