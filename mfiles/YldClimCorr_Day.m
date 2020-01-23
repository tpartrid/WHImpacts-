
clear; close all
path(path,'/ihome/tpartrid/Matlab/Toolboxes/')

load /adata/tpartrid/TBagger_ProcData/Input/County/TB_InputFeatures.mat
load /adata/tpartrid/NASS/County/CGY_BUA_AllCNTY.mat
load /adata/tpartrid/NASS/County/CGAH_CNTY_1930_2015.mat
load /adata/tpartrid/TBagger_ProcData/Input/County/IncludedCounties.mat
%     S = shaperead('/adata/tpartrid/NASS/Shapefiles/County/cb_2016_us_county_20m.shp');

%Exclude Yield from coutnies with less than 10000 acres harvested 
cols = IncludedCnty;     %Limit Yield to  districts with >30000 acres harvested on average 
CGY_AllCnty(:,cols~=1) = NaN; 
    
    Txu_mx_d1 = NaN(1,size(CGY_AllCnty,2));
    Txu_mx_d2 = NaN(1,size(CGY_AllCnty,2));
    Txx_mx_d1 = NaN(1,size(CGY_AllCnty,2));
    Txx_mx_d2 = NaN(1,size(CGY_AllCnty,2));
    Tnu_mx_d1 = NaN(1,size(CGY_AllCnty,2));
    Tnu_mx_d2 = NaN(1,size(CGY_AllCnty,2));
    Tnn_mx_d1 = NaN(1,size(CGY_AllCnty,2));
    Tnn_mx_d2 = NaN(1,size(CGY_AllCnty,2));
    Pr_mx_d1 = NaN(1,size(CGY_AllCnty,2));
    Pr_mx_d2 = NaN(1,size(CGY_AllCnty,2));
    Prcorr_out = NaN(1,size(CGY_AllCnty,2));
    Txucorr_out = NaN(1,size(CGY_AllCnty,2)); Txxcorr_out = NaN(1,size(CGY_AllCnty,2));
    Tnucorr_out = NaN(1,size(CGY_AllCnty,2)); Tnncorr_out = NaN(1,size(CGY_AllCnty,2));

%     PlntDate = SoyPlntDate;
    for c = 1:size(CGY_AllCnty,2)
        c
       
        Yield = CGY_AllCnty(:,c); % SoyY_AgDist(1:86,c);%average yield over district
        yr = (1930:2015)'; 
        %ID missing data in Yield and Temperature - some counties have
        %yield data but no temp values from CDL mask
        md_Y = ~isnan(Yield); md_T = ~isnan(Tx_gs(:,c));
        Y = NaN(size(Yield,1),1);
        if sum(md_Y) > 30 && sum(md_T)>0
            %Detrend Yield 
            [curve, goodness, output] = fit(yr(md_Y),Yield(md_Y),'smoothingspline','SmoothingParam',.001); %'smoothingspline','SmoothingParam',0.05
            Y(md_Y) = output.residuals;
            tic
            
            for y = 1930:2015
                idx = t(:,2)==y;
                GDD_temp =GDD(idx,c);
                KDD_temp =KDD(idx,c);
                tn_temp =Tn_gs(idx,c);
                tx_temp =Tx_gs(idx,c);
                pr_temp =Pr_gs(idx,c);
                
                
                gdd(1:gsl(c)+1,y-1930+1) = GDD_temp(~isnan(GDD_temp)); 
                kdd(1:gsl(c)+1,y-1930+1) = KDD_temp(~isnan(KDD_temp));
                tn(1:gsl(c)+1,y-1930+1) = tn_temp(~isnan(tn_temp)); 
                tx(1:gsl(c)+1,y-1930+1) = tx_temp(~isnan(tx_temp)); 
                pr(1:gsl(c)+1,y-1930+1) = pr_temp(~isnan(pr_temp)); 
            end
            

            for d1 = 1:size(gdd,1)
                d1;
                for d2 = d1:size(gdd,1)
                    d2;
                    gddin = nansum(gdd(d1:d2,:),1)';
                    kddin = nansum(kdd(d1:d2,:),1)';
                            
                    prin = nansum(pr(d1:d2,:),1)';
                    tnin = nanmean(tn(d1:d2,:),1)';
                    txin = nanmean(tx(d1:d2,:),1)';
                    txxin = nanmax(tx(d1:d2,:),[],1)';
                    tnnin = nanmin(tn(d1:d2,:),[],1)';

                     [r,p] = corr(Y,prin,'rows','complete');
                     Prcorr(d1,d2) = r;
                     Prcorr_p(d1,d2) = p;
                     
                     [r,p] = corr(Y,tnin,'rows','complete');
                     Tnucorr(d1,d2) = r;
                     Tnucorr_p(d1,d2) = p;
                     
                     [r, p] = corr(Y,txin,'rows','complete');
                     Txucorr(d1,d2) = r;
                     Txucorr_p(d1,d2) = p;
                     
                     [r,p] = corr(Y,tnnin,'rows','complete');
                     Tnncorr(d1,d2) = r;
                     Tnncorr_p(d1,d2) = p;
                     
                     [r, p] = corr(Y,txxin,'rows','complete');
                     Txxcorr(d1,d2) = r;
                     Txxcorr_p(d1,d2) = p;
                     
                     [r,p] = corr(Y,gddin,'rows','complete');
                     GddCorr(d1,d2) = r;
                     GddCorr_p(d1,d2) = p;
                     
                     [r,p] = corr(Y,kddin,'rows','complete');
                     KDDCorr(d1,d2) = r;
                     KDDCorr_p(d1,d2) = p;

                end
            end
        
            clear gdd kdd pr tx tn
            Txucorr(Txucorr_p >= 0.01) = NaN;
            Tnucorr(Tnucorr_p >= 0.01) = NaN;
            Txxcorr(Txxcorr_p >= 0.01) = NaN;
            Tnncorr(Tnncorr_p >= 0.01) = NaN;
            Prcorr(Prcorr_p >= 0.01) = NaN;
            GddCorr(GddCorr_p >= 0.01) = NaN;
            KDDCorr(KDDCorr_p >= 0.01) = NaN;
            
            hf1 = figure; %Show correlation table as image 
            dk1 = 1:gsl(c)+1;
            dk2 = 1:gsl(c)+1;
            imagesc(dk2,dk1,Txxcorr)
            hold on 
            cl = colorbar; 
            caxis([-0.6 0.6])
                cl.Label.String = '\rho';
            set(gca,'FontSize',20)
            xlabel('Last Day')
            ylabel('First Day')
            title(['Daily correlation: Tx & Corn, ' num2str(c)])

%             %Make animation of correlation plot for each district to see if
%             mov(c) = getframe(hf1);
%             export_fig Txcorr.tif -append -nocrop
%             im2gif('Txcorr.tif','Txcorr_Corn.gif','-delay',1/8,'-nocrop')
%             close 
%             
% 
%             hf2 = figure; %Show correlation table as image 
%             dk1 = 1:gsl(c)+1;
%             dk2 = 1:gsl(c)+1;
%             imagesc(dk2,dk1,Prcorr)
%             hold on 
%             cl = colorbar; 
%             caxis([-0.6 0.6])
%                 cl.Label.String = '\rho';
%             set(gca,'FontSize',20)
%             xlabel('Last Day')
%             ylabel('First Day')
%             title(['Daily correlation: Pr & Corn, ' num2str(c)])
% 
%             %Make animation of correlation plot for each district to see if
%             mov2(c) = getframe(hf2);
%             export_fig Prcorr.tif -append -nocrop
%             im2gif('Prcorr.tif','Prcorr_Corn.gif','-delay',1/8,'-nocrop')
%             close 

            %Find min or max correlation- i.e. most important period for yield 
            %max - finds only positive correlation i.e. beneficial 
            %min - find negative correlation i.e. harmful
            %max(abs(val - doesn't assume direciton... 
            %Anyway to account for drought directly ??? 
           
            %Prcp 
            [mn1, wkidx] = max(Prcorr);
            [~, Prwk2] = max(mn1);
            Pr_wk2(c) = Prwk2;
            Pr_wk1(c) = wkidx(Pr_wk2(c));
            Prcorr_out(c) = Prcorr(Pr_wk1(c),Pr_wk2(c));

            %Txu
            [mn1, wkidx] = min(Txucorr); %max(abs ??
            [~, Txuwk2] = min(mn1);
            Txu_wk2(c) = Txuwk2;
            Txu_wk1(c) = wkidx(Txu_wk2(c)); 
            Txucorr_out(c) = Txucorr(Txu_wk1(c),Txu_wk2(c));

            %Tnu
            [mn1, wkidx] = min(Tnucorr);
            [~, Tnuwk2] = min(mn1);
            Tnu_wk2(c) = Tnuwk2;
            Tnu_wk1(c) = wkidx(Tnu_wk2(c));
            Tnucorr_out(c) = Tnucorr(Tnu_wk1(c),Tnu_wk2(c));
            
            %Txx
            [mn1, wkidx] = min(Txxcorr); 
            [~, Txxwk2] = min(mn1);
            Txx_wk2(c) = Txxwk2;
            Txx_wk1(c) = wkidx(Txx_wk2(c)); 
            Txxcorr_out(c) = Txxcorr(Txx_wk1(c),Txx_wk2(c));

            %Tnn
            [mn1, wkidx] = min(Tnncorr); % Assume that Tnn is also harmful 
            [~, Tnnwk2] = min(mn1);
            Tnn_wk2(c) = Tnnwk2;
            Tnn_wk1(c) = wkidx(Tnn_wk2(c));
            Tnncorr_out(c) = Tnncorr(Tnn_wk1(c),Tnn_wk2(c));

            %GDD
            [mx1, wkidx] = min(GddCorr); %max(abs ??
            [~, Txuwk2] = min(mx1);
            GDD_d2(c) = Txuwk2;
            GDD_d1(c) = wkidx(GDD_d2(c)); 
            GDDCorr_out(c) = GddCorr(GDD_d1(c),GDD_d2(c));
            %KDD
            [mx1, wkidx] = min(KDDCorr); %max(abs ??
            [~, Txuwk2] = min(mx1);
            KDD_d2(c) = Txuwk2;
            KDD_d1(c) = wkidx(KDD_d2(c)); 
            KDDCorr_out(c) = KDDCorr(KDD_d1(c),KDD_d2(c));

%         %Temporary 
%             %Prcp 
%             [mx1, wkidx] = max(Prcorr);
%             [~, Prwk2] = max(mx1);
%             Pr_mx_wk2(c) = Prwk2;
%             Pr_mx_wk1(c) = wkidx(Pr_mx_wk2(c));
%             Prcorr_out(c) = max(mx1);
%             %PRCP min
%             [mn1, wkidx] = min(Prcorr);
%             [~, Prwk2] = min(mn1);
%             Pr_mn_wk2(c) = Prwk2;
%             Pr_mn_wk1(c) = wkidx(Pr_mn_wk2(c));
%             Prcorr_mn_out(c) = min(mn1);
% 
%             %Txu
%             [mx1, wkidx] = max(Txucorr); %max(abs ??
%             [~, Txuwk2] = max(mx1);
%             Txu_mx_wk2(c) = Txuwk2;
%             Txu_mx_wk1(c) = wkidx(Txu_mx_wk2(c)); 
%             Txucorr_out(c) = max(mx1);
%             %Tx Min
%             [mn1, wkidx] = min(Txucorr); %max(abs ??
%             [~, Txuwk2] = min(mn1);
%             Txu_mn_wk2(c) = Txuwk2;
%             Txu_mn_wk1(c) = wkidx(Txu_mn_wk2(c)); 
%             Txucorr_mn_out(c) = min(mn1);
% 
%             %Tnu
%             [mx1, wkidx] = max(Tnucorr);
%             [~, Tnuwk2] = max(mx1);
%             Tnu_mx_wk2(c) = Tnuwk2;
%             Tnu_mx_wk1(c) = wkidx(Tnu_mx_wk2(c));
%             Tnucorr_out(c) = max(mx1);
%             %Tnu Min
%             [mn1, wkidx] = min(Tnucorr);
%             [~, Tnuwk2] = min(mn1);
%             Tnu_mn_wk2(c) = Tnuwk2;
%             Tnu_mn_wk1(c) = wkidx(Tnu_mn_wk2(c));
%             Tnucorr_mn_out(c) = min(mn1);
%             
%             %Txx
%             [mx1, wkidx] = max(Txxcorr); 
%             [~, Txxwk2] = max(mx1);
%             Txx_mx_wk2(c) = Txxwk2;
%             Txx_mx_wk1(c) = wkidx(Txx_mx_wk2(c)); 
%             Txxcorr_out(c) = max(mx1);
%             %Txx Min
%             [mn1, wkidx] = min(Txxcorr); 
%             [~, Txxwk2] = min(mn1);
%             Txx_mn_wk2(c) = Txxwk2;
%             Txx_mn_wk1(c) = wkidx(Txx_mn_wk2(c)); 
%             Txxcorr_mn_out(c) = min(mn1);
% 
%             %Tnn
%             [mx1, wkidx] = max(Tnncorr); % Assume that Tnn is also harmful 
%             [~, Tnnwk2] = max(mx1);
%             Tnn_mx_wk2(c) = Tnnwk2;
%             Tnn_mx_wk1(c) = wkidx(Tnn_mx_wk2(c));
%             Tnncorr_out(c) = max(mx1);
%             %Tnn Min
%             [mn1, wkidx] = min(Tnncorr); % Assume that Tnn is also harmful 
%             [~, Tnnwk2] = min(mn1);
%             Tnn_mn_wk2(c) = Tnnwk2;
%             Tnn_mn_wk1(c) = wkidx(Tnn_mn_wk2(c));
%             Tnncorr_mn_out(c) = min(mn1);
%             clear Prcorr Txucorr Tnucorr Txxcorr Tnncorr mn1 mn2 wkidx Prwk2 Tnnwk2 Txxwk2 Tnuwk2 Txxwk2 
        end
        clear Txucorr Tnucorr Prcorr GddCorr KDDCorr Txxcorr Tnncorr Txucorr_p Tnucorr_p Txxcorr_p Tnncorr_p Prcorr_p KDDCorr_p GddCorr_p
        toc
    end

%     CritPeriod_mx = cat(1, Txu_mx_wk1, Txu_mx_wk2, Txx_mx_wk1, Txx_mx_wk2, Tnu_mx_wk1, Tnu_mx_wk2, Tnn_mx_wk1, Tnn_mx_wk2, Pr_mx_wk1, Pr_mx_wk2); 
%     CritPeriod_mn = cat(1, Txu_mn_wk1, Txu_mn_wk2, Txx_mn_wk1, Txx_mn_wk2, Tnu_mn_wk1, Tnu_mn_wk2, Tnn_mn_wk1, Tnn_mn_wk2, Pr_mn_wk1, Pr_mn_wk2); 
CritPeriod = cat(1, Txu_wk1, Txu_wk2, Tnu_wk1, Tnu_wk2, Pr_wk1, Pr_wk2, GDD_d1, GDD_d2,KDD_d1,KDD_d2, Txx_wk1, Txx_wk2, Tnn_wk1, Tnn_wk2); 
%     CritPeriodVars = {'TxuWk1';'TxuWk2';'TxxWk1';'TxxWk2';'TnuWk1';'TnuWk2';'TnnWk1',;'TnnWk2';'PrWk1';'PrWk2'}; 
CritPeriodVars = {'Txud1';'Txud2';'Tnud1';'Tnud2';'PrWk1';'PrWk2';'GDDd1';'GDDd2';'KDDday1';'KDDday2';'Txxd1';'Txxd2';'Tnnd1';'Tnnd2'}; 
CritPeriod(CritPeriod == 0) = NaN;

% save /adata/tpartrid/TBagger_ProcData/Input/County/YieldClimCorr_Corn_day.mat Txu_wk1, Txu_wk2, Tnu_wk1, Tnu_wk2, Pr_wk1, Pr_wk2, GDD_d1, GDD_d2,KDD_d1,KDD_d2
save /adata/tpartrid/TBagger_ProcData/Input/County/YieldClimCorr_Corn_day2_corrected.mat CritPeriod CritPeriodVars Txucorr_out Tnucorr_out Prcorr_out GDDCorr_out KDDCorr_out
toc