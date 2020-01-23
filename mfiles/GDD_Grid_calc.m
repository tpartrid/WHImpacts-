%GDD and KDD gridded calc - claculates GDD and KDD for each grid in the
%Livneh daataset 

%Livneh Files for interpolation  - Trim while loading to save memory 
% load /adata/tpartrid/Livneh_1915_2015/latlon_Livneh.mat; %Livneh grid 
% rows = lat >= 36 & lat <= 48; %Add 1° buffer on either sdie 
% cols = lon >= -103 & lon <= -81;
% load /adata/tpartrid/CultivatedArea.mat
%  load /adata/tpartrid/Livneh_1915_2015/Scenarios/latlon_LivnehTrim.mat
load /adata/tpartrid/Livneh_1915_2015/Scenarios/Tn_Orig.mat
% z = t(:,2)>=1930 & t(:,2) <= 2015; 
% Tn_dly_Delta(rows~=1,:,:) = []; Tn_dly_Delta(:,cols~=1,:) = []; Tn_dly_Delta(:,:,z~=1) = [];

load /adata/tpartrid/Livneh_1915_2015/Scenarios/Tx_Orig.mat
z = t(:,2)>=1930 & t(:,2) <= 2011; 
% Tx_dly_Delta(rows~=1,:,:) = []; Tx_dly_Delta(:,cols~=1,:) = []; Tx_dly_Delta(:,:,z~=1) = [];
% t(z~=1,:) = [];
% lat(rows~=1) = [];
% lon(cols~=1) = [];
% save /adata/tpartrid/Livneh_1915_2015/Tx_Livneh_CROP.mat Tx_dly_Delta t lat lon -v7.3 
% save /adata/tpartrid/Livneh_1915_2015/Tn_Livneh_CROP.mat Tn_dly_Delta t lat lon -v7.3 

% load /adata/tpartrid/Livneh_1915_2015/Tx_Livneh_CROP.mat 
% load /adata/tpartrid/Livneh_1915_2015/Tn_Livneh_CROP.mat 
%Growing Degree Days, Killing Degree Days, and Average GS Tx, Tn, Pr 
Tlow = 9; %Constant in GDD calculation -  9°C (based on Butler and Huybers) 
Thigh = 29; %Constant in GDD calculation - 29


%Growing Degree Days 
%Tmax_d = Tx_dly_QMAP;
%Tmax_d(Tmax_d<=Tlow) = Tlow;
%Tmax_d(Tmax_d>=Thigh) = Thigh;
%Tn
%Tmin_d = Tn_dly_QMAP;
%Tmin_d(Tmin_d<=Tlow) = Tlow;
%Tmin_d(Tmin_d>=Thigh) = Thigh;
GDD =(Tx_dly_Orig+Tn_dly_Orig)/2 - 8;
GDD(GDD<=0) = 0; 
clear Tmin_d Tn_dly_Orig Tmax_d

%Killing Degree Days 
%Tmax_d = Tx_dly_QMAP;
Tmax_d = Tx_dly_Orig - 33;
Tmax_d(Tmax_d <= 0) = 0;
KDD = Tmax_d;

% %For Looping it ... 
% for i = 1:size(Tx_dly_Delta,1)
%     i
%     for j = 1:size(Tx_dly_Delta,2)
%         j
%         for d = 1:size(Tx_dly_Delta,3)
%             d
%             %Determine Tx Daily 
%             if Tlow < Tx_dly_Delta(i,j,d) && Tx_dly_Delta(i,j,d) < Thigh %TxDaily = orig value if between Tlow and THigh
%                 TxDaily = Tx_dly_Delta(i,j,d)
%             elseif Tx_dly_Delta(i,j,d) <= Tlow %TxDaily = Tlow if less than Tlow
%                 TxDaily = Tlow
%             elseif Tx_dly_Delta(i,j,d) >= Thigh %TxDaily = Thigh is greater than Thigh
%                 TxDaily = Thigh
%             end
%             %Determine Tn Daily - Same as for TxDaily
%             if Tlow < Tn_dly_Delta(i,j,d) && Tn_dly_Delta(i,j,d) < Thigh
%                 TnDaily = Tn_dly_Delta(i,j,d);
%             elseif Tn_dly_Delta(i,j,d) <= Tlow
%                 TnDaily = Tow;
%             elseif Tn_dly_Delta(i,j,d) >= Thigh
%                 TnDaily = Thigh;
%             end
%             %Calculate Daily GDD
%             GDDtemp = ((TxDaily +TnDaily)/2) - Tlow;
%             %Calculate KDD : Tmax - Thigh
%             if Txtemp(i,j,d) > Thigh
%                 KDDtemp = Txtemp(i,j,d) - Thigh; %KDD = Daily Tx - Thigh unless Daily Tx is less than Thigh
%             else  %Txtemp(d) <= Thigh
%                 KDDtemp = 0;
%             end
% 
%             GDD(i,j,d) = GDDtemp; %Store yearly values in full array 
%             KDD(i,j,d) = KDDtemp;
%     %         GDD_tot((i,j,d-1930+1,j) = nansum(GDDtemp); %Calculate cumulative GDD and KDD values 
%     %         KDD_tot((i,j,d-1930+1,j) = nansum(KDDtemp);
%         end
%         clear KDDtemp GDDtemp
%     end
% end

save /adata/tpartrid/Livneh_1915_2015/Scenarios/GDDKDD_grid_Orig.mat KDD GDD t -v7.3 
