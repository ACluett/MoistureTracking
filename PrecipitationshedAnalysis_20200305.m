%%%%%% SET LOCATION %%%%%%

clear

cd '/Users/allisoncluett/Dropbox/GitHub/MoistureTracking'

% Set location HERE and file names will change accordingly throughout
location = "THULE";

if location == 'KANG'
    
    % Set site latitude (deg N) and longitude (deg E)
    target_lat = 61.229;
    target_lon = -48.097;
    
    % Set GHCN station identifier
    station = 'GLE00147241';
    
else if location  == 'SUAQ'
        
    % Set site latitude (deg N) and longitude (deg E)
    target_lat = 69.21666667
    target_lon = -51.05
    
    % Set GHCN station identifier
    station = 'GLE00146961';

    else if location == 'THULE'
            
      % Set site latitude (deg N) and longitude (deg E)
      target_lat = 77.467;
      target_lon = -69.228;
      
      % Set GHCN station identifier
      station = 'GLW00017605';
      
        end
    end
end


%% Read in neccessary files


% Load tracked evaporation data
E_track_file = strcat(location,"_evaporation_tracking_data/E_track_2layers_sigma_continental2014-1979.mat");
load(E_track_file);

% Load mask file of sink region
region_file = strcat(location,"_region_grid.mat");
load(region_file);

% Load color scheme for plotting maps
load("20200121_RevisedColormap.mat");

% Load spatial data
load("Grid.mat");

% Load custom region masks
load("Custom_region_masks.mat");

% Load other reanalysis and GHCN climate data
load("ERA-Int_Temp_Formatted.mat");
GHCN_All = readtable("GHCN.csv");

% Calculate area of sink region (area grid is in m2)
region_area = sum(sum(Region .* full_area_grid));

% Crop off first and last years from tracked evaporation
E_track_per_year_per_month_raw = E_track_per_year_per_month;
E_track_per_year_per_month = E_track_per_year_per_month(2:35, :, :, :);

% Crop off first and last years track reanalysis precipitation
P_per_year_per_month_raw = P_per_year_per_month;
P_per_year_per_month = P_per_year_per_month(2:35, :, :, :);


%% Compare amount of precipitation between station data, reanalysis precip, and tracked evap

% Pull out GHCN data for site out of table which can contain multiple sites
SiteTable = GHCN_All(strcmp(GHCN_All.STATION, station), :);

% Make a vector for each variable (check that columns are ordered the same)
GHCN_PRCP = SiteTable{:,8};
GHCN_TEMP = SiteTable{:,9};
GHCN_MONTH = SiteTable{:,7};
GHCN_YEAR = SiteTable{:,6};

% Count number of observations in monthly GHCN precip and temp records

for i = 1:12 
PRCP_N(i) = sum(~isnan(GHCN_PRCP(GHCN_MONTH==i)));
end

for i = 1:12 
TEMP_N(i) = sum(~isnan(GHCN_TEMP(GHCN_MONTH==i)));
end

% Calculate  monthly averages and standard deviation for GHCN temp
GHCN_Temp_MonthAvg = zeros(1,12);
GHCN_Temp_MonthStd = zeros(1,12);

for Midx = 1:12
    GHCN_Temp_MonthAvg(Midx) = mean(GHCN_TEMP(GHCN_MONTH==Midx), 'omitnan');
    GHCN_Temp_MonthStd(Midx) = std(GHCN_TEMP(GHCN_MONTH==Midx), 'omitnan');
end

% Calculate monthly totals for reanalysis precip and tracked evaporation

for Yidx = 1:34
    for Midx = 1:12
        
    % 1000 multiplier is convert from m to mm (unit of GHCN data)
    Etrack_month_sum(Yidx, Midx) = squeeze(sum(sum(E_track_per_year_per_month(Yidx, Midx, :, :)))) / region_area * 1000;
    P_slice = squeeze(P_per_year_per_month(Yidx,Midx,:,:));
    reanalysis_P(Yidx, Midx) = sum(sum(P_slice .* Region))/region_area * 1000;
    
    end
end

% Calculate annual sums for GHCN data
% Make a vector of years for which the GHCN data covers
GHCN_unique_years = unique(SiteTable.YEAR);

for Yidx = 1:length(GHCN_unique_years)
    
    GHCN_unique_year_data = SiteTable(SiteTable.YEAR == GHCN_unique_years(Yidx),:);
    
    % Check if each year has 12 months of data
    % if it does add them up and save 'em, if not NaN
    
    if sum(~isnan(GHCN_unique_year_data.PRCP)) == 12
        GHCN_annual_sum(Yidx) = sum(GHCN_unique_year_data.PRCP, 1);
    else
        GHCN_annual_sum(Yidx) = NaN;
        
    end
end

% Count sumber of observations in GHCN data
GHCN_annual_n = sum(~isnan(GHCN_annual_sum));

% Calculate annual sums for reanalysis precip and Etrack
Etrack_annual_sum = sum(Etrack_month_sum, 2);
reanalysisP_annual_sum = sum(reanalysis_P, 2);


% Plot boxplots to compare GHCN precip, reanalysis precip, and tracked evap at annual scale

figure
subplot(1,3,1)
boxplot(Etrack_annual_sum)
title('Tracked Evap')
ylabel('Annual precip (mm)')
ylim([0 2200])
o = findobj(gcf,'tag','Outliers');
set(o,'Marker', 'o','MarkerSize',3)

subplot(1,3,2)
boxplot(reanalysisP_annual_sum)
title('Reanalysis Precip')
ylabel('Annual precip (mm)')
ylim([0 2200])
o2 = findobj(gcf,'tag','Outliers');
set(o2,'Marker', 'o','MarkerSize',3)

subplot(1,3,3)
boxplot(GHCN_annual_sum)
title('GHCN Precip')
ylabel('Annual precip (mm)')
ylim([0 2200])
o3 = findobj(gcf,'tag','Outliers');
set(o3,'Marker', 'o','MarkerSize',3)

% Uncomment the following lines to save in specific directory
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 15 5]) 
% plota1=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Precip_Validation/',datestr(now, 'yyyy-mm-dd'),location,'AnnualPrecip_Boxplots.eps');
% saveas(gcf,plota1)


% Run a 2 sample t-test to identify months with significantly different mean precip between datasets

    for Midx = 1:12
    
    month_reanalysisP = reanalysis_P(:, Midx);
    month_Etrack = Etrack_month_sum(:, Midx);
    month_GHCN = GHCN_PRCP(GHCN_MONTH == Midx);

    GHCN_Etrack_ttest(Midx) = ttest2(month_GHCN, month_Etrack)
    GHCN_reanalysisP_ttest(Midx) = ttest2(month_GHCN, month_reanalysisP)
    Etrack_reanalysisP_ttest(Midx) = ttest2(month_Etrack, month_reanalysisP)

    end

% Run a 2 sample t-test to compare annual sums
GHCN_Etrackttest_annual = ttest2(GHCN_annual_sum, Etrack_annual_sum)
    
% Calculate and plot residual between reanalysis P and tracked evaporation
reanalysis_residual = reanalysis_P - Etrack_month_sum;
figure
errorbar(mean(reanalysis_residual), std(reanalysis_residual))
title("Precipitation Residual")
xlabel("Month")
ylabel("Precipitation-Tracked Evaporation (mm)")
xlim([0 13]);
ylim([-50 50]);

% Uncomment the following lines to save in specific directory
%plotc=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Precip_Validation/',datestr(now, 'yyyy-mm-dd'),location,'Precip_Residuals.eps');
%saveas(gcf,plotc)


%% Plot annual average map for absolute values

% Set latitude and longitude limits for map
latlim= [0 90];
lonlim= [-180 180];

% Make a geographic reference matrix
geoR = makerefmat('RasterSize', [108 240], ...
'Latlim', [min(latitude) max(latitude)] , 'Lonlim', [min(longitude) max(longitude)]);

fig=figure;
axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,...
'Frame','on','Grid', 'on','MeridianLabel','on','ParallelLabel','on')
h = meshm(flip(squeeze(sum(mean(E_track_per_year_per_month)))./region_area  .* 1000),geoR);
set(gca,'ColorScale','linear')
colormap(gca, mymap);
c=colorbar;

% Scale colorbar for location
if location == 'KANG'
    caxis([0 8]);  
    
else if location  == 'SUAQ'
    caxis([0 10]);  
    
    else if location == 'THULE'
       caxis([0 1]);     
            
        end
    end
end

c.Location = 'southoutside';
set(get(c,'XLabel'),'String','Contribution to Mean Annual Site Precip (mm)')
plotm(target_lat , target_lon,'kp', 'MarkerSize', 8, 'MarkerFaceColor','k')
load coastlines
plotm(coastlat,coastlon, 'k')


axis off

%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])  
%plot0=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Annual_Avg/',datestr(now, 'yyyy-mm-dd'),location,'Annual_Avg.png');
%saveas(gcf,plot0)
%To export illustrator editable fig: need to manually change file name
%print -painters -depsc annualavg_Thule

%% Plot annual average map for fractions

% Calculate array of fractional contribution from each grid cell
nYears = 34;
nMonths = 12;
nlat = length(latitude);
nlon = length(longitude);

% Set up array to store fractions of annual tracked evap
fraction_array = zeros(34,12,108,240);

% Set up array to store fractions of tracked evap for each timeslice
slice_fraction_array = zeros(34,12,108,240);

% Calculate fractional contribution from each grid cell in each timeslice
nYears = 34;
nMonths = 12;
nLat = 108;
nLon = 240;

for Yidx = 1:nYears
    
    year_sum = sum(sum(sum(E_track_per_year_per_month(Yidx,:,:,:))));
    
    for Midx = 1:nMonths
       
        slice = squeeze(squeeze(E_track_per_year_per_month(Yidx, Midx,:,:)));
        slice_sum = sum(sum(slice));
        
        for Latidx = 1:nLat
            for Lonidx = 1 :nLon
                fraction_array(Yidx, Midx, Latidx, Lonidx) = E_track_per_year_per_month(Yidx, Midx, Latidx, Lonidx)/year_sum;
                slice_fraction_array(Yidx, Midx, Latidx, Lonidx) = E_track_per_year_per_month(Yidx, Midx, Latidx, Lonidx)/slice_sum;
            end
        end
        
    end
end

latlim= [0 90];
lonlim= [-180 30];

geoR = makerefmat('RasterSize', [108 240], ...
'Latlim', [min(latitude) max(latitude)] , 'Lonlim', [min(longitude) max(longitude)]);

fig=figure;
axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,...
'Frame','on','Grid', 'on','MeridianLabel','on','ParallelLabel','on')
h = meshm(flip(squeeze(sum(mean(fraction_array))).*100),geoR);
set(gca,'ColorScale','linear')
colormap(gca, mymap);
c=colorbar;
caxis([0 0.45]);     
            
c.Location = 'southoutside';
set(get(c,'XLabel'),'String','Percent Contribution to Mean Annual Site Precip (mm)')
plotm(target_lat , target_lon,'kp', 'MarkerSize', 8, 'MarkerFaceColor','k')
load coastlines
plotm(coastlat,coastlon, 'k')

axis off

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3])  
% plot0=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Annual_Avg/',datestr(now, 'yyyy-mm-dd'),location,'Annual_Avg.png');
% saveas(gcf,plot0)
%To export illustrator editable fig: need to manually change file name
%print -painters -depsc annualfractions_KANG

%% Plot fractional grid cell contributions for all months in one figure

latlim= [0 90];
lonlim= [-180 30];

geoR = makerefmat('RasterSize', [108 240], ...
'Latlim', [min(latitude) max(latitude)] , 'Lonlim', [min(longitude) max(longitude)]);

figure;

% Set up spacing to tile all months in one figure
set(gcf,'position',[103 438 2354 907],'PaperPositionMode','auto','BackingStore','off','PaperOrientation','landscape');
cnt=0;
nr=3;
nc=4;
mn=0;
mx=0.5;

for jj=1:12
    cnt=cnt+1;
    s(jj)=subplot(nr,nc,cnt);
end

cnt=0;

for jj=1:12
    cnt=cnt+1;
    subplot(s(jj))
    if rem(jj,nc)==0
        s(jj).Position=[1-(1/nc) 1-(ceil(jj/nc)/nr) 1/nc 1/nr];
    else
        s(jj).Position=[(rem(jj,nc)/nc)-(1/nc) 1-(ceil(jj/nc)/nr) 1/nc 1/nr];
    end
    
        axesm('robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,...
              'Frame','on','Grid', 'on','MeridianLabel','off','ParallelLabel','off')
        h = meshm(flip(squeeze(mean(slice_fraction_array(:,jj,:,:)))).* 100 ,geoR);
        set(gca,'ColorScale','linear')
        colormap(gca, mymap);
        c=colorbar;
        caxis([0 0.45]);
        colorbar('off')

        load coastlines
        plotm(coastlat,coastlon, 'k')
        axis off
        
      drawnow;

end

% Uncomment to print illustrator-editable figure
%print -painters -depsc monthlyfractions_SUAQv2



%% Calculate regional sums and percents for each time slice

region_sums = zeros(34, 12, 16);
region_percents = zeros(34, 12, 16);
nYears = 34;
nMonths = 12;
nRegions = 16;

for Yidx = 1:nYears
    
    for Midx = 1:nMonths
       
        slice = squeeze(squeeze(E_track_per_year_per_month(Yidx, Midx,:,:)));
        
        region_sums(Yidx, Midx, 1) = sum(sum(mask_Other .* slice))/region_area;
       
        region_sums(Yidx, Midx, 2) = sum(sum(mask_NAmN .* slice)) /region_area;
        region_sums(Yidx, Midx, 3) = sum(sum(mask_NAMnC .* slice))/region_area;
        region_sums(Yidx, Midx, 4) = sum(sum(mask_NAmC .* slice))/region_area;
        region_sums(Yidx, Midx, 5) = sum(sum(mask_NAmS .* slice))/region_area;
        
        region_sums(Yidx, Midx, 6) = sum(sum(mask_ArcticOcean .* slice))/region_area;
        region_sums(Yidx, Midx, 7) = sum(sum(mask_GINSea .* slice))/region_area;
        
        region_sums(Yidx, Midx, 8) = sum(sum(mask_CanadianArctic .* slice))/region_area;
        region_sums(Yidx, Midx, 9) = sum(sum(mask_Baffin .* slice))/region_area;
        region_sums(Yidx, Midx, 10) = sum(sum(mask_LabSea .* slice))/region_area;
        
        region_sums(Yidx, Midx, 11) = sum(sum(mask_AtlanticN .* slice))/region_area;
        region_sums(Yidx, Midx, 12) = sum(sum(mask_AtlanticC .* slice))/region_area;
        region_sums(Yidx, Midx, 13) = sum(sum(mask_AtlanticS .* slice))/region_area;
        
        region_sums(Yidx, Midx, 14) = sum(sum(mask_PacificN .* slice))/region_area;
        region_sums(Yidx, Midx, 15) = sum(sum(mask_PacificC .* slice))/region_area;
        region_sums(Yidx, Midx, 16) = sum(sum(mask_PacificS .* slice))/region_area;
        
        region_percents(Yidx, Midx, 1) = region_sums(Yidx, Midx, 1) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 2) = region_sums(Yidx, Midx, 2) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 3) = region_sums(Yidx, Midx, 3) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 4) = region_sums(Yidx, Midx, 4) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 5) = region_sums(Yidx, Midx, 5) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 6) = region_sums(Yidx, Midx, 6) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 7) = region_sums(Yidx, Midx, 7) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 8) = region_sums(Yidx, Midx, 8) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 9) = region_sums(Yidx, Midx, 9) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 10) = region_sums(Yidx, Midx, 10) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 11) = region_sums(Yidx, Midx, 11) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 12) = region_sums(Yidx, Midx, 12) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 13) = region_sums(Yidx, Midx, 13) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 14) = region_sums(Yidx, Midx, 14) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 15) = region_sums(Yidx, Midx, 15) / sum(region_sums(Yidx, Midx, :));
        region_percents(Yidx, Midx, 16) = region_sums(Yidx, Midx, 16) / sum(region_sums(Yidx, Midx, :));
        
    end

end


%% Plot monthly average regional sum bar chart

monthly_region_avgs = zeros(12, 16);

for Midx = 1:nMonths
    for Ridx = 1:nRegions
    monthly_region_avgs(Midx, Ridx) = mean(region_sums(:,Midx,Ridx))* 1000;
    end
end

figure('Position', [10 10 500 500]);

H = bar(monthly_region_avgs, 'stacked');

% Define custom color scheme
set(H(1),'facecolor', '#323030')

set(H(2),'facecolor', '#5E3B42')
set(H(3),'facecolor', '#78565A')
set(H(4),'facecolor', '#A78E8D')
set(H(5),'facecolor', '#C8B7B4')

set(H(10),'facecolor', '#F3EBD2')
set(H(9),'facecolor', '#E0CC8F')
set(H(8),'facecolor', '#CCAC52')
set(H(7),'facecolor', '#948140')
set(H(6),'facecolor', '#7E6924')

set(H(11),'facecolor', '#578484')
set(H(12),'facecolor', '#95B4B1')
set(H(13),'facecolor', '#AABABA')

set(H(14),'facecolor', '#6B7E9D')
set(H(15),'facecolor', '#8998B1')
set(H(16),'facecolor', '#B7C1CC')

ax = gca; 
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.LineWidth = 1;
box on
ax.XColor = 'k';
ax.YColor = 'k';
lgd.TextColor = 'k';
legend('Location','southoutside');
lgd.NumColumns = 2;
ldg.EdgeColor = 'k';
ylim([0 100]);
xlabel('Month')
ylabel('Summed E-track (mm)')

if location == 'KANG'
    str = 'A) Kangilinnguit';
    else if location == 'SUAQ'
            str = 'B) Kangerlussuaq';
        else if location == 'THULE'
                str = 'C) Thule';
            end
        end
end

annotation('textbox', [0.13 .955 0 0], 'String', str, 'FitBoxToText', true, 'LineStyle', 'none', 'FontSize', 14);

legend('Other', 'North American: North', 'North America: North Central', 'North America: Central', 'North America: South', ...
     'Arctic Ocean', 'GIN Seas', 'Canadian Arctic', 'Baffin Bay', 'Labrador Sea', 'Atlantic: North', ...
         'Atlantic: Central', 'Atlantic: South', 'Pacific: North', 'Pacific: Central', 'Pacific: South')

% plot1=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Sector_Division/',datestr(now, 'yyyy-mm-dd'),location,'SectorSums.svg');
% saveas(gcf,plot1)

%% Plot monthly average regional percent bar chart

monthly_region_percent_avgs = zeros(12, 16);

for Midx = 1:nMonths
    for Ridx = 1:nRegions
    monthly_region_percent_avgs(Midx, Ridx) = mean(region_percents(:,Midx,Ridx));
    end
end

figure('Position', [10 10 500 500]);

J = bar(monthly_region_percent_avgs, 'stacked');

set(J(1),'facecolor', '#323030')

set(J(2),'facecolor', '#5E3B42')
set(J(3),'facecolor', '#78565A')
set(J(4),'facecolor', '#A78E8D')
set(J(5),'facecolor', '#C8B7B4')

set(J(10),'facecolor', '#F3EBD2')
set(J(9),'facecolor', '#E0CC8F')
set(J(8),'facecolor', '#CCAC52')
set(J(7),'facecolor', '#948140')
set(J(6),'facecolor', '#7E6924')

set(J(11),'facecolor', '#578484')
set(J(12),'facecolor', '#95B4B1')
set(J(13),'facecolor', '#AABABA')

set(J(14),'facecolor', '#6B7E9D')
set(J(15),'facecolor', '#8998B1')
set(J(16),'facecolor', '#B7C1CC')

ax = gca; 
ax.FontSize = 12;
ax.TickDir = 'out';
ax.TickLength = [0.01 0.01];
ax.LineWidth = 1;
box on
ax.XColor = 'k';
ax.YColor = 'k';
lgd.TextColor = 'k';
legend('Location','southoutside');
lgd.NumColumns = 2;
ldg.EdgeColor = 'k';

xlabel('Month')
ylabel('Percent E-track (m)')

if location == 'KANG'
    str = 'A) Kangilinnguit';
    else if location == 'SUAQ'
            str = 'B) Kangerlussuaq';
        else if location == 'THULE'
                str = 'C) Thule';
            end
        end
end


annotation('textbox', [0.13 .955 0 0], 'String', str, 'FitBoxToText', true, 'LineStyle', 'none', 'FontSize', 14);

legend('Other', 'North American: North', 'North America: North Central', 'North America: Central', 'North America: South', ...
     'Arctic Ocean', 'GIN Seas', 'Canadian Arctic', 'Baffin Bay', 'Labrador Sea', 'Atlantic: North', ...
         'Atlantic: Central', 'Atlantic: South', 'Pacific: North', 'Pacific: Central', 'Pacific: South')

plot2=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Sector_Division/',datestr(now, 'yyyy-mm-dd'),location,'SectorPercents.svg');
saveas(gcf,plot2)

%Save matrix of monthly region percents
% data1 = strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Sector_Division/',datestr(now, 'yyyy-mm-dd'),location,'SectorPercents.mat')
% save(data1, 'monthly_region_percent_avgs');

%% Make matrix of distances of grid cell center to target location

Dist_Mat = zeros(108, 240);

for i = 1:length(latitude)
    for j = 1:length(longitude)
                
        Dist_Mat(i,j) = deg2km(distance([target_lat, target_lon], [latitude(i) + 0.75, longitude(j) - 0.75]));
        
    end
end

%% Calculate transport distance and average source latitude for each timeslice
% Then calculate monthly averages and stdevs

sea_mask = ~lsm;

slice_avg_dist = [];
slice_avg_dist_land = [];
slice_avg_dist_sea = [];
slice_fr_land = [];
slice_fr_sea = [];
month_avg_dist = [];
month_avg_dist_land = [];
month_avg_dist_sea = [];
slice_dist_mat = [];

latitude_mat = repmat(latitude,1, 240);

for Midx = 1:12
    
   for Yidx = 1: 34
          
        slice = squeeze(E_track_per_year_per_month(Yidx,Midx,:,:));
        
        slice_sum = sum(sum(slice));
        slice_fr = slice ./ slice_sum;
        slice_dist_weights = slice_fr .* Dist_Mat;
        slice_avg_dist(Yidx) = sum(sum(slice_dist_weights));
        slice_dist_mat(Yidx, Midx) = sum(sum(slice_dist_weights));
        
        slice_avg_lat(Yidx) = sum(sum(slice_fr .* latitude_mat));
        
        landslice = slice .* lsm;
        slice_landsum = sum(sum(landslice));
        slice_amount_weights_land = landslice ./ slice_landsum;
        slice_dist_weights_land = slice_amount_weights_land .* Dist_Mat;
        slice_lat_weights_land = slice_amount_weights_land .* latitude_mat;
        slice_avg_dist_land(Yidx) = sum(sum(slice_dist_weights_land));
        slice_avg_lat_land(Yidx) = sum(sum(slice_lat_weights_land));
        
        seaslice = slice .* sea_mask;
        slice_seasum = sum(sum(seaslice));
        slice_amount_weights_sea = seaslice ./ slice_seasum;
        slice_dist_weights_sea = slice_amount_weights_sea .* Dist_Mat;
        slice_lat_weights_sea = slice_amount_weights_sea .* latitude_mat;
        slice_avg_dist_sea(Yidx) = sum(sum(slice_dist_weights_sea));
        slice_avg_lat_sea(Yidx) = sum(sum(slice_lat_weights_sea));
        
        slice_fr_land(Yidx) = slice_landsum/slice_sum;
        slice_fr_sea(Yidx) = slice_seasum/slice_sum;

   end
   
   month_avg_dist(Midx) = mean(slice_avg_dist);
   month_std_dist(Midx) = std(slice_avg_dist);
   month_avg_dist_land(Midx) = mean(slice_avg_dist_land);
   month_avg_dist_sea(Midx) = mean(slice_avg_dist_sea);
   
   month_avg_lat(Midx) = mean(slice_avg_lat);
   month_std_lat(Midx) = std(slice_avg_lat);
   month_avg_lat_land(Midx) = mean(slice_avg_lat_land);
   month_avg_lat_sea(Midx) = mean(slice_avg_lat_sea);
   
   month_fr_land(Midx) = mean(slice_fr_land);
   month_fr_sea(Midx) = mean(slice_fr_sea);
   
   
end

figure

subplot(3,1,1)
hold on
plot(month_fr_land, '--r')
plot(month_fr_sea, '--b')
hold off
xlim([0 13])
ylim([0 1])
ylabel('Fraction Moisture')

subplot(3,1,2)

hold on

errorbar(month_avg_dist, month_std_dist, '-k')
plot(month_avg_dist_sea,  '-b')
plot(month_avg_dist_land, '-r')
ylim([2000 6000])
xlim([0 13])

ylabel('Weighted Distance (km)')

lgd.TextColor = 'k';
legend('Location','northeast');
lgd.NumColumns = 2;
ldg.EdgeColor = 'k';

subplot(3,1,3)
hold on
errorbar(month_avg_lat, month_std_lat, '-k')
plot(month_avg_lat_sea, 'b')
plot(month_avg_lat_land, 'r')
hold off
xlim([0 13])
ylim([30 60])
xlabel('Month')
ylabel('Latitude (°N)')

plot6=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Transport_Dist/',datestr(now, 'yyyy-mm-dd'),location,'AvgTransportDistStd.eps');
saveas(gcf,plot6)

% To export for looking at interannual variability/correlations
% dist_file = strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/Transport_Dist/',datestr(now, 'yyyy-mm-dd'),location,'TotalTransportDist.csv');
% writematrix(slice_dist_mat, dist_file)


%% Calculate monthly cooling between source and sink

slice_avg_T2m_source = zeros(34, 12);
slice_avg_T2m_source_land = zeros(34, 12);
slice_avg_T2m_source_sea = zeros(34, 12);
slice_avg_T2m_sink = zeros(34, 12);
slice_avg_T2m_instrumental_sink = zeros(34,12);
slice_avg_D2m_sink = zeros(34, 12);
slice_avg_D2m_instrumental_sink = zeros(34,12);
slice_avg_T2m_diff = zeros(34, 12);
slice_dist_mat = zeros(34, 12);
slice_avg_T2m_fixedsource = zeros(34, 12);
slice_avg_T2m_fixedseason = zeros(34, 12);
slice_avg_SST_source_sea = zeros(34,12);
slice_avg_T2m_diff_fixedsource = zeros(34,12);
slice_avg_T2m_source_fixedseason = zeros(34,12);
slice_avg_T2m_diff_fixedseason = zeros(34,12);

month_avg_T2m_sink = zeros(1,12);
month_std_T2m_sink = zeros(1,12);
month_avg_D2m_sink = zeros(1,12);

month_avg_T2m_instrumental_sink = zeros(1,12);
month_avg_D2m_instrumental_sink = zeros(1,12);

month_avg_T2m_source = zeros(1,12);
month_avg_T2m_source_land = zeros(1,12);
month_avg_T2m_source_sea = zeros(1,12);
month_avg_T2m_source_fixedsource = zeros(1,12);
month_avg_T2m_source_fixedseason = zeros(1,12);
month_avg_T2m_diff = zeros(1,12);
month_avg_D2m_diff = zeros(1,12);
month_std_T2m_diff = zeros(1,12);
month_std_D2m_diff = zeros(1,12);
month_avg_SST = zeros(1,12);


for Midx = 1:12
    
   for Yidx = 1: 34
          
        Etrack_slice = squeeze(E_track_per_year_per_month(Yidx,Midx,:,:));
        T2m_slice = squeeze(T2m(Yidx, Midx,:,:));
        D2m_slice = squeeze(D2m(Yidx, Midx,:,:));
        SST_slice = squeeze(SST(Yidx, Midx,:,:));
        Area_slice_sink = full_area_grid .* Region;
        slice_avg_T2m_sink(Yidx, Midx) = sum((sum(Area_slice_sink .* T2m_slice)))./ region_area;
        slice_avg_D2m_sink(Yidx, Midx) = sum((sum(Area_slice_sink .* D2m_slice)))./ region_area;
        T2m_slice_fixedseason = squeeze(squeeze(mean(mean(T2m(:,:,:,:)))));
        
        Etrack_slice_sum = sum(sum(Etrack_slice));
        Etrack_slice_fr = Etrack_slice ./ Etrack_slice_sum;
        T2m_slice_weights_source = Etrack_slice_fr .* T2m_slice;
        slice_avg_T2m_source(Yidx, Midx) = sum(sum(T2m_slice_weights_source));
        
        slice_avg_T2m_diff(Yidx, Midx) = slice_avg_T2m_source(Yidx, Midx) - slice_avg_T2m_sink(Yidx, Midx);
        slice_avg_D2m_diff(Yidx, Midx) = slice_avg_T2m_source(Yidx, Midx) - slice_avg_D2m_sink(Yidx, Midx);
        
        Etrack_landslice = Etrack_slice .* double(lsm);
        Etrack_slice_sum_land = sum(sum(Etrack_landslice));
        Etrack_slice_fr_land = Etrack_landslice ./ Etrack_slice_sum_land;
        T2m_slice_weights_source_land = Etrack_slice_fr_land .* T2m_slice;
        slice_avg_T2m_source_land(Yidx, Midx) = sum(sum(T2m_slice_weights_source_land));

        SST_mask = ~isnan(SST_slice);
        
        Etrack_seaslice = Etrack_slice .* SST_mask;
        Etrack_slice_sum_sea = sum(sum(Etrack_seaslice));
        Etrack_slice_fr_sea = Etrack_seaslice ./ Etrack_slice_sum_sea;
        T2m_slice_weights_source_sea = Etrack_slice_fr_sea .* T2m_slice;
        slice_avg_T2m_source_sea(Yidx, Midx) = sum(sum(T2m_slice_weights_source_sea));
        
        SST_slice_weights_source_sea = Etrack_slice_fr_sea .* SST_slice;
        slice_avg_SST_source_sea(Yidx, Midx) = sum(sum(SST_slice_weights_source_sea, 'omitnan'));
        
        Etrack_slice_fixedsource = squeeze(squeeze(mean(mean(E_track_per_year_per_month(:,:,:,:)))));
        Etrack_slice_sum_fixedsource = sum(sum(Etrack_slice_fixedsource));
        Etrack_slice_fr_fixedsource = Etrack_slice_fixedsource ./ Etrack_slice_sum_fixedsource;
        T2m_slice_weights_fixedsource = Etrack_slice_fr_fixedsource .* T2m_slice;
       
        slice_avg_T2m_source_fixedsource(Yidx, Midx) = sum(sum(T2m_slice_weights_fixedsource));
        slice_avg_T2m_diff_fixedsource(Yidx, Midx) = slice_avg_T2m_source_fixedsource(Yidx, Midx) - slice_avg_T2m_sink(Yidx, Midx);
        
        slice_avg_T2m_source_fixedseason(Yidx, Midx) = sum(sum(Etrack_slice_fr .* T2m_slice_fixedseason));
        slice_avg_T2m_diff_fixedseason(Yidx, Midx) = slice_avg_T2m_source_fixedseason(Yidx, Midx) - slice_avg_T2m_sink(Yidx, Midx);
        
   end
   
   % Monthly averages for plots
   
   month_avg_T2m_sink(Midx) = mean(slice_avg_T2m_sink(:,Midx)) -273.15;
   month_std_T2m_sink(Midx) = std(slice_avg_T2m_sink(:, Midx));
   
   month_avg_D2m_sink(Midx) = mean(slice_avg_D2m_sink(:,Midx))-273.15;
   
   month_avg_T2m_instrumental_sink(Midx) = mean(slice_avg_T2m_instrumental_sink(:,Midx))-273.15;
   month_avg_D2m_instrumental_sink(Midx) = mean(slice_avg_D2m_instrumental_sink(:,Midx))-273.15;
   
   month_avg_T2m_source(Midx) = mean(slice_avg_T2m_source(:,Midx))-273.15;
   month_avg_T2m_source_land(Midx) = mean(slice_avg_T2m_source_land(:,Midx))-273.15;
   month_avg_T2m_source_sea(Midx) = mean(slice_avg_T2m_source_sea(:,Midx))-273.15;
   month_avg_T2m_source_fixedsource(Midx) = mean(slice_avg_T2m_source_fixedsource(:,Midx)) - 273.15;
   month_avg_T2m_source_fixedseason(Midx) = mean(slice_avg_T2m_source_fixedseason(:,Midx)) -273.15;
   
   month_avg_T2m_diff(Midx) = mean(slice_avg_T2m_diff(:, Midx));
   month_avg_D2m_diff(Midx) = mean(slice_avg_D2m_diff(:, Midx));
   month_std_T2m_diff(Midx) = std(slice_avg_T2m_diff(:, Midx));
   month_std_D2m_diff(Midx) = std(slice_avg_D2m_diff(:, Midx));
   
   month_avg_SST(Midx) = mean(slice_avg_SST_source_sea(:, Midx))-273.15;
   
end


%% Simple Rayleigh Distillation Model using GHCN monthly mean temps as Tf

% Define three scenarios TEMPERATURES IN KELVIN
month_avg_T2m_sinkK = month_avg_T2m_sink + 273.15;
month_avg_T2m_sourceK = month_avg_T2m_source + 273.15;
month_avg_T2m_source_fixedsourceK = month_avg_T2m_source_fixedsource + 273.15;
month_avg_T2m_source_fixedseasonK = month_avg_T2m_source_fixedseason + 273.15;

Tf = month_avg_T2m_sinkK;

% Scenario 1: Local seasonal temp change, source temp constant
T0_1 = repmat(mean(month_avg_T2m_sourceK), 12, 1);

% Scenario 2: Local seasonal temp change, fixed source seasonal temp change
T0_2 = month_avg_T2m_source_fixedsourceK;

% Scenario 3: Local seasonal temp change, source dist and seasonal temp change
T0_3 = month_avg_T2m_sourceK;

% % % Scenario 4:Local seasonal temp change, moving source with no seasonal temp change
T0_4 = month_avg_T2m_source_fixedseasonK;

%%%%%% Run all 4 scenarios for site %%%%%%%%

% Pick a starting isotopic composition
 
% Constant scenario

d2Hv0a = [-120 -120 -120 -120 -120 -120 -120 -120 -120 -120 -120 -120];

% Vary based on oceanic source latitude (d2Hv0b)

if location == 'KANG'
    
    d2Hv0b = [-118.82 -118.75 -113.96 -109.17 -106.00 -96.61 -94.97 -99.54 -110.80 -112.96 -114.30 -119.05];
    
    else if location == 'SUAQ'
            
                d2Hv0b = [-119.37 -118.76 -114.38 -112.70 -111.29 -98.27 -96.96 -105.22 -115.17 -116.63 -117.87 -119.98];
                
        else if location == 'THULE'
                   
                d2Hv0b = [-119.76 -119.38 -115.28 -113.01 -113.81 -102.87 -100.36 -111.49 -118.70 -120.04 -119.36 -120.46];                  
            end
            
        end
        
end


for i = 1:12

% Scenario 1

fm_1(i) = vapor(Tf(i))/vapor(T0_1(i));
d2Hvf_1(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_1(i)));
d2Hrf_1(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_1(i))) + E2Hpv(Tf(i));


% Scenario 2

fm_2(i) = vapor(Tf(i))/vapor(T0_2(i));
d2Hvf_2(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_2(i)));
d2Hrf_2(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_2(i))) + E2Hpv(Tf(i));


% Scenario 3

fm_3(i) = vapor(Tf(i))/vapor(T0_3(i));

d2Hvf_3a(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_3(i)));
d2Hrf_3a(i) = d2Hv0a(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_3(i))) + E2Hpv(Tf(i));

d2Hvf_3b(i) = d2Hv0b(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_3(i)));
d2Hrf_3b(i) = d2Hv0b(i) + E2Hpv(283.15) .* log(vapor(Tf(i))/vapor(T0_3(i))) + E2Hpv(Tf(i));




end

%% Plot Rayleigh distillation results

figure

subplot(3,1,1)

hold on

plot(Tf-273.15,  '-vk', 'MarkerFaceColor', 'k')

plot(T0_1-273.15, '--^k', 'MarkerFaceColor', 'r')
plot(T0_2-273.15, '-^k', 'MarkerFaceColor', 'k')
plot(T0_3-273.15, '-^k', 'MarkerFaceColor', 'b')


hold off

ylim([-50 30])
xlim([0 13])
ylabel('Temp (°C)')

subplot(3,1,2)
hold on

plot(fm_1, 'r')
plot(fm_2, 'k')
plot(fm_3, 'b')

hold off

xlim([0 13])
ylim([0 1])
ylabel('fm')

subplot(3,1,3)

hold on
plot(d2Hrf_1, 'r')
plot(d2Hrf_2, 'k')
plot(d2Hrf_3a, '-b')
plot(d2Hrf_3b, '--b')
ylabel('d2H Precip (? VSMOW)')

hold off

ylim([-400 50])
xlim([0 13])

lgd.TextColor = 'k';
legend('Location','southeast');
lgd.NumColumns = 2;
ldg.EdgeColor = 'k';
legend('Scenario 1: Fixed Source Annual T0', 'Scenario 2: Seasonal Fixed Source T0', 'Scenario 3a: Moving Source Seasonal T0 Constant Vapor', 'Scenario 3b: latitudinal seasonal oceanic vapor change')

hold off

% plot9=strcat('/Users/allisoncluett/Dropbox/UB/Research/Dissertation/Ch.2 Moisture Variability/Precipitationsheds/precip_Validation/',datestr(now, 'yyyy-mm-dd'),location,'RayleighSimulation_TfGHCN_DewpointCrxn.eps');
% saveas(gcf,plot9)



%% Define functions for use in Rayleigh model

% Calculate saturation vapor pressure as a function of temperature (in K)

function V = vapor(T)
V = 0.0002.*(T-273.15).^3 + 0.0111.*(T-273.15).^2 + 0.321.*(T-273.15) + 4.8;
end

% Temperature dependent fractionation between vapor and precipitation
% (liquid or ice)

function E2H = E2Hpv(T)
if T > 273.15;
    E2H = (24.844.*(10^6./T.^2)) + (-76.248 .* (10.^3./T)) + 52.612; % (Majoube, 1971)
else T <= 273.15;
    E2H = (24.844.*(10^6./T.^2)) + (-76.248 .* (10.^3./T)) + 71.912;
    end
end



