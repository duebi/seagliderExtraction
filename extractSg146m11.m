% extractSg146m11.m
%
% script to extract data from 2015 SeaGlider mission sg146m11
% 
% Benedetto Barone - Oct 2015

mission = 'sg146_m11';
d_range = [2 321]; 
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
clear upth
load([sgpath '/oxy_cal'])
load([sgpath '/hplc_cal'])

%load pcpn_cal
cd(sgpath)
nd = length(d_range(1):d_range(2)); % total number of dives
%{
% 1. Collect data from all dives:
[UP,DWN,~] = sg_cat(d_range,sgpath,'userVarNames',{'vmtime','press','sigmath0','salin', ...
    'tempc','oxygen','optode_oxygen','optode_temp','optode_dphase_oxygen','optode_dphase','wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'});

% 2. Rename WETLabs variables
DWN.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb660','cdom','chl2'};
UP.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb660','cdom','chl2'};

% 3. Clean backscattering data based on 3 std on a 20m grid (sets outliers = NaN)
for i = d_range(1):d_range(2)
    ind_d = DWN.divenum==i & ~isnan(DWN.bb470);
    [bin_d,bin,ind_good_d] = binning(DWN.bb470(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb470(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb470(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb470);
    [bin_u,bin,ind_good_u] = binning(UP.bb470(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb470(ind_u); temp_u(~ind_good_u) = NaN; UP.bb470(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb660);
    [bin_d,bin,ind_good_d] = binning(DWN.bb660(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb660(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb660(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb660);
    [bin_u,bin,ind_good_u] = binning(UP.bb660(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb660(ind_u); temp_u(~ind_good_u) = NaN; UP.bb660(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb700);
    [bin_d,bin,ind_good_d] = binning(DWN.bb700(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb700(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb700(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb700);
    [bin_u,bin,ind_good_u] = binning(UP.bb700(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb700(ind_u); temp_u(~ind_good_u) = NaN; UP.bb700(ind_u) = temp_u;
    %{
    plot(UP.bb470(ind_u),UP.vmdepth(ind_u),'k.',DWN.bb470(ind_d),DWN.vmdepth(ind_d),'r.')
    hold on, plot(bin_u,bin,'k',bin_d,bin,'r'), hold off
    set(gca,'ydir','rev','ylim',[0 200])
    pause
    %}
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
end

% 4. Compute particle backscattering coefficients (this could be moved into sg_divecalc.m)
chi_p = 1.1;
ldwn = height(DWN); lup = height(UP); 
betasw470d = NaN(ldwn,1); betasw660d = NaN(ldwn,1); betasw700d = NaN(ldwn,1);
betasw470u = NaN(lup,1); betasw660u = NaN(lup,1); betasw700u = NaN(lup,1);
    % compute scatter at 117 deg due to seawater
for i = 1:ldwn % downcast 
    [betasw470d(i),~,~]= betasw_ZHH2009(470,DWN.tempc(i),117,DWN.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw660d(i),~,~]= betasw_ZHH2009(660,DWN.tempc(i),117,DWN.salin(i));
    [betasw700d(i),~,~]= betasw_ZHH2009(700,DWN.tempc(i),117,DWN.salin(i));
end
for i = 1:lup % upcast
    [betasw470u(i),~,~]= betasw_ZHH2009(470,UP.tempc(i),117,UP.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw660u(i),~,~]= betasw_ZHH2009(660,UP.tempc(i),117,UP.salin(i));
    [betasw700u(i),~,~]= betasw_ZHH2009(700,UP.tempc(i),117,UP.salin(i));
end
    % calculation for particle backscattering coefficients
DWN.bbp470 = 2*pi*chi_p*(DWN.bb470-betasw470d);
DWN.bbp660 = 2*pi*chi_p*(DWN.bb660-betasw660d);
DWN.bbp700 = 2*pi*chi_p*(DWN.bb700-betasw700d);
UP.bbp470 = 2*pi*chi_p*(UP.bb470-betasw470u);
UP.bbp660 = 2*pi*chi_p*(UP.bb660-betasw660u);
UP.bbp700 = 2*pi*chi_p*(UP.bb700-betasw700u);
clear ldwn lup sgpath
clear chi_p betasw470d betasw470u betasw660d betasw660u betasw700d betasw700u

% 5. Correct optode oxygen for salinity and pressure (because dphase is missing in this mission)
    % Aandera salinity correction coefficients (sensor independent)
aa_sb = [-6.24097E-3 -6.93498E-3 -6.90358E-3 -4.29155E-3];
aa_sc = -3.11680E-7;
    % Salinity correction
scal_t = log((298.15 - DWN.optode_temp)./(273.15 + DWN.optode_temp));
DWN.optode_oxygen =  DWN.optode_oxygen.*exp(DWN.salin.*(aa_sb(1)+aa_sb(2)*scal_t+aa_sb(3)*scal_t.^2+aa_sb(4)*scal_t.^3)+aa_sc.*DWN.salin.^2);
scal_t = log((298.15 - UP.optode_temp)./(273.15 + UP.optode_temp));
UP.optode_oxygen =  UP.optode_oxygen.*exp(UP.salin.*(aa_sb(1)+aa_sb(2)*scal_t+aa_sb(3)*scal_t.^2+aa_sb(4)*scal_t.^3)+aa_sc.*UP.salin.^2);
    % Pressure correction
DWN.optode_oxygen = DWN.optode_oxygen.*(1+DWN.press.*0.032./1000);
UP.optode_oxygen = UP.optode_oxygen.*(1+UP.press.*0.032./1000);

% 6. Save concatenated dive file
save([mission '_cat'],'UP','DWN');
%}


% DOWNCAST
% 7. Put data on a regular grid
datafile = [mission(1:5) mission(7)  mission(8:9) 'data'];
depth = 2:2:1000; ld = length(depth);
[UG,DG] = sg_grid(['./' mission '_cat'],depth,'gridVar','vmdepth','diveRange',d_range(1):d_range(2),'outVars',{'lon','lat','daten','salin','tempc','sigmath0','oxygen','optode_dphase_oxygen','optode_oxygen','chl1','chl2','bbp470','bbp660','bbp700','cdom'});
    % Build new table for downcast
    % Time & Position
sgd = table(depth',reshape(DG.daten,ld,nd)); sgd.Properties.VariableNames = ({'depth','date'});
sgd.date = sgd.date-10/24; %to transform in HST time
sgd.hour = sgd.date - fix(sgd.date);
sgd.lon = reshape(DG.lon,ld,nd);
sgd.lat = reshape(DG.lat,ld,nd);
    % Water column measurements
sgd.s = reshape(DG.salin,ld,nd);
sgd.t = reshape(DG.tempc,ld,nd);
sgd.sig = reshape(DG.sigmath0,ld,nd);
sgd.o = reshape(DG.oxygen,ld,nd).*((sgd.sig+1000)/1000); % From umol kg-1 to umol L-1
sgd.opt = reshape(DG.optode_oxygen,ld,nd); % The dphase is missing for this cruise, so I had to use the data that have been computed internally
sgd.chl1 = reshape(DG.chl1,ld,nd);
sgd.chl2 = reshape(DG.chl2,ld,nd);
sgd.bbp470 = reshape(DG.bbp470,ld,nd);
sgd.bbp660 = reshape(DG.bbp660,ld,nd);
sgd.bbp700 = reshape(DG.bbp700,ld,nd);
sgd.cdom = reshape(DG.cdom,ld,nd);
sgd.Properties.VariableUnits = {'m','Matlab date (HST)','Decimal day (HST)','Degrees E','Degrees N'...
    'g kg-1 (Absolute Salinity)', 'degrees Celsius', 'kg m-3','umol L-1','umol L-1','rel. units','rel. units' ...
    'm-1','m-1','m-1', 'rel. units'};
    % Dive information
dived.lon = nanmean(sgd.lon);
dived.lat = nanmean(sgd.lat);
dived.dive = d_range(1):d_range(2);
dived.date = nanmean(sgd.date);
dived.hour = dived.date - fix(dived.date);
[dived.u10,~,~,~,~,~] = ow_2obs(dived.lat,dived.lon,dived.date+10/24); % extract wind data
dived.u10 = dived.u10';
dived.u10(isnan(dived.u10)) = interp1(dived.date(~isnan(dived.u10)),dived.u10(~isnan(dived.u10)),dived.date(isnan(dived.u10))); % interpolate missing data
[dived.slp,~,~,~] = nr_2obs(dived.lat,dived.lon,dived.date+10/24,'pres','surface_gauss'); % Sea level pressure in Pascals
dived.slp = dived.slp';
dived.units = {'lon: Degrees E';'lat: Degrees N';'dive: #';'date: Matlab date (HST)';'hour: Decimal day (HST)';'u10: wind at 10 m (m s-1)';'slp: sea level pressure (Pa)'};

% 8. Oxygen calibration on Winkler measurements
    % build comparison array
o2_comp = oxy_cal;
o2_comp.dist = NaN(height(oxy_cal),1);
o2_comp.lag = NaN(height(oxy_cal),1);
o2_comp.zdist = NaN(height(oxy_cal),1);
o2_comp.optode = NaN(height(oxy_cal),1);
o2_comp.sbe43 = NaN(height(oxy_cal),1);
for i = 1:height(oxy_cal)
    [o2_comp.lag(i),ind_t] = min(abs(dived.date-o2_comp.date(i)));
    o2_comp.dist(i) = vdist(dived.lat(ind_t),dived.lon(ind_t),o2_comp.lat(i),-o2_comp.lon(i));
    [o2_comp.zdist(i),ind_z] = min(abs(sgd.depth-o2_comp.depth(i)));
    o2_comp.optode(i) = sgd.opt(ind_z,ind_t);
    o2_comp.sbe43(i) = sgd.o(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10;
ind_correction = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10 & o2_comp.depth<30;
    % Compute offset in the upper 100 m and correct oxygen from sensors
optode_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction));
sbe43_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction));
sgd.opt = sgd.opt + optode_offset;
sgd.o = sgd.o + sbe43_offset;
    % Plot comparison results
%{
subplot(2,2,1)
scatter(o2_comp.oxygen(ind_good),o2_comp.optode(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('Optode oxygen (umol L-1)')
subplot(2,2,2)
scatter(o2_comp.oxygen(ind_good),o2_comp.sbe43(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('SBE43 oxygen (umol L-1)')
subplot(2,2,3)
hist(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - optode (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
subplot(2,2,4)
hist(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - SBE43 (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
%}

% 9. Fluorescence-Chlorophyll calibration from HPLC pigment measurements
% build comparison array
hplc_comp = hplc_cal(:,{'DATE','DEPTH','LAT','LON','TOTCHLA'});
hplc_comp.TPP = hplc_cal.TOTCHLA + hplc_cal.CHLC + hplc_cal.MV_CHLB + hplc_cal.FUCO + hplc_cal.HEX + hplc_cal.BUT;
hplc_comp.dist = NaN(height(hplc_cal),1);
hplc_comp.lag = NaN(height(hplc_cal),1);
hplc_comp.zdist = NaN(height(hplc_cal),1);
hplc_comp.chl1 = NaN(height(hplc_cal),1);
hplc_comp.chl2 = NaN(height(hplc_cal),1);
for i = 1:height(hplc_cal)
    [hplc_comp.lag(i),ind_t] = min(abs(dived.date-hplc_comp.DATE(i)));
    hplc_comp.dist(i) = vdist(dived.lat(ind_t),dived.lon(ind_t),hplc_comp.LAT(i),-hplc_comp.LON(i));
    [hplc_comp.zdist(i),ind_z] = min(abs(sgd.depth-hplc_comp.DEPTH(i)));
    hplc_comp.chl1(i) = sgd.chl1(ind_z,ind_t);
    hplc_comp.chl2(i) = sgd.chl2(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = hplc_comp.dist<20000 & hplc_comp.lag<0.2 & hplc_comp.zdist<10;
    % Compute offset in the upper 100 m and correct oxygen from sensors
% Regression fluorescence vs. chla+chlb+chlc+psc
[param_c1,stat_c1] = robustfit(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good));
[slope_psc1,intercept_psc1,r_psc1,sslope_psc1,sintercept_psc1]  = lsqfitgm(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good));
[param_c2,stat_c2] = robustfit(hplc_comp.chl2(ind_good),hplc_comp.TPP(ind_good));
opts = optimset('MaxIter',20000,'MaxFunEval',200000);
beta3 = lsqcurvefit('gompertzfun',[0.89 -0.31 -3000 -0.08],hplc_comp.DEPTH,hplc_comp.TOTCHLA./hplc_comp.TPP,[],[],opts);
r_sig3 = corrcoef(gompertzfun(beta3,hplc_comp.DEPTH),hplc_comp.TOTCHLA./hplc_comp.TPP);
% Calibrate fluorescence on chlorophyll a
sgd.chl1 = (param_c1(2)*sgd.chl1 + param_c1(1)).*(gompertzfun(beta3,repmat(sgd.depth,1,nd)));
sgd.chl2 = (param_c2(2)*sgd.chl2 + param_c2(1)).*(gompertzfun(beta3,repmat(sgd.depth,1,nd)));
sgd.Properties.VariableUnits({'chl1','chl2'}) = {'ng chl-a L-1','ng chl-a L-1'};
    % Plot comparison results
%{
subplot(2,2,1)
plot(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good),'ko')
hold on, plot([nanmin(hplc_comp.chl1) nanmax(hplc_comp.chl1)],param_c1(2)*[nanmin(hplc_comp.chl1) nanmax(hplc_comp.chl1)]+param_c1(1),'k--'), hold off
xlabel('Seaglider Fluorescence 1'), ylabel('Total Photosynthetic Pigments(ng/L)'), set(gca,'Fontsize',18)
lg = legend('data','robust fit'), set(lg,'box','off','Location','NorthWest','Fontsize',18)
subplot(2,2,3)
plot(hplc_comp.chl2(ind_good),hplc_comp.TPP(ind_good),'ko')
hold on, plot([nanmin(hplc_comp.chl2) nanmax(hplc_comp.chl2)],param_c2(2)*[nanmin(hplc_comp.chl2) nanmax(hplc_comp.chl2)]+param_c2(1),'k--'), hold off
xlabel('Seaglider Fluorescence 2'), ylabel('Total Photosynthetic Pigments(ng/L)'), set(gca,'Fontsize',18)
subplot(2,2,[2 4])
plot(hplc_comp.TOTCHLA./hplc_comp.TPP,hplc_comp.DEPTH,'ko'),hold on, plot(gompertzfun(beta3,0:200),0:1:200,'k--'), hold off, set(gca,'ydir','rev','Fontsize',18)
xlabel('Chlorophyll a / Total Photosynthetic Pigments(ng/ng)'), ylabel('Depth (m)')
%}
    
% Compute characteristics on isopycnal levels (average sigma at depth bins)
    % compute average sigma at depth bins
sig_grid = nanmean(sgd.sig,2);
sig_grid = sort(sig_grid);
    % initialize table of isopycnal data
isod = sgd;
isod.sig = sig_grid;
isod.date = NaN*isod.date; isod.hour = NaN*isod.hour; isod.lat = NaN*isod.lat; isod.lon = NaN*isod.lon;
isod.s = NaN*isod.s; isod.t = NaN*isod.t; isod.o = NaN*isod.o; isod.opt = NaN*isod.opt;
isod.chl1 = NaN*isod.chl1; isod.chl2 = NaN*isod.chl2; isod.cdom = NaN*isod.cdom;
isod.bbp470 = NaN*isod.bbp470; isod.bbp660 = NaN*isod.bbp660; isod.bbp700 = NaN*isod.bbp700;
    % interpolate at isopycnal values
for i = 1:nd
    ind_nan = isnan(sgd.sig(:,i));
    if sum(~ind_nan) ~= 0
        isod.t(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.t(~ind_nan,i),sig_grid);
        isod.s(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.s(~ind_nan,i),sig_grid);
        isod.date(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.date(~ind_nan,i),sig_grid);
        isod.hour(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.hour(~ind_nan,i),sig_grid);
        isod.o(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.o(~ind_nan,i),sig_grid);
        isod.opt(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.opt(~ind_nan,i),sig_grid);
        isod.chl1(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.chl1(~ind_nan,i),sig_grid);
        isod.chl2(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.chl2(~ind_nan,i),sig_grid);
        isod.cdom(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.cdom(~ind_nan,i),sig_grid);
        isod.bbp470(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.bbp470(~ind_nan,i),sig_grid);
        isod.bbp660(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.bbp660(~ind_nan,i),sig_grid);
        isod.bbp700(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.bbp700(~ind_nan,i),sig_grid);
    end
    clear ind_nan
end
    
% Compute surface oxygen flux
[Fd,Fp,Fc,~] = fas_L13(sgd.opt(2,:)/1000,dived.u10,sgd.s(2,:)*35/35.16504,sgd.t(2,:),dived.slp*0.00000986923267,'O2');
dived.Fs = (Fd+Fp+Fc)*86400; % from s-1 to d-1
dived.units{end+1} = 'Fs: surface oxygen flux (mol m-2 d-1)';
clear Fd Fp Fc

% InterpolateSLA on seaglider position and time
[dived.sla,~,~,~] = av_2obs(dived.lat,dived.lon,dived.date+10/24,'sla','linear');
dived.sla = dived.sla*100;
dived.units{end+1} = 'sla: sea level anomaly from AVISO (cm)';

% Find mixed layer depth and density at the base of the mixed layer (0.03 kg m-3 difference from sigma at 10 m depth)
mld = sgd.sig - repmat(sgd.sig(5,:),height(sgd),1) - 0.03;
mld(mld<0) = NaN; mld(1,:) = NaN;
[sig003,ind003] = nanmin(mld);
mld003 = sgd.depth(ind003);
mld003sig = sig003 + sgd.sig(2,:) + 0.03;
dived.mld003 = mld003';
dived.mld003sig = mld003sig;
dived.units{end+1} = 'mld003: mixed layer depth with drho 0.03 kg m-3 from 10m (m)';
dived.units{end+1} = 'mld003sig: potential density anomaly at MLD (kg m-3)';
clear mld ind003 sig003 mld003 mld003sig  

% UPCAST
    % Build new table for upcast
    % Time & Position
sgu = table(depth',reshape(UG.daten,ld,nd)); sgu.Properties.VariableNames = ({'depth','date'});
sgu.date = sgu.date-10/24; %to transform in HST time
sgu.hour = sgu.date - fix(sgu.date);
sgu.lon = reshape(UG.lon,ld,nd);
sgu.lat = reshape(UG.lat,ld,nd);
    % Water column measurements
sgu.s = reshape(UG.salin,ld,nd);
sgu.t = reshape(UG.tempc,ld,nd);
sgu.sig = reshape(UG.sigmath0,ld,nd);
sgu.o = reshape(UG.oxygen,ld,nd).*((sgu.sig+1000)/1000); % From umol kg-1 to umol L-1
sgu.opt = reshape(UG.optode_oxygen,ld,nd); % The dphase is missing for this cruise, so I had to use the data that have been computed internally
sgu.chl1 = reshape(UG.chl1,ld,nd);
sgu.chl2 = reshape(UG.chl2,ld,nd);
sgu.bbp470 = reshape(UG.bbp470,ld,nd);
sgu.bbp660 = reshape(UG.bbp660,ld,nd);
sgu.bbp700 = reshape(UG.bbp700,ld,nd);
sgu.cdom = reshape(UG.cdom,ld,nd);
sgu.Properties.VariableUnits = {'m','Matlab date (HST)','Decimal day (HST)','Degrees E','Degrees N'...
    'g kg-1 (Absolute Salinity)', 'degrees Celsius', 'kg m-3','umol L-1','umol L-1','rel. units','rel. units' ...
    'm-1','m-1','m-1', 'rel. units'};
    % Dive information
diveu.lon = nanmean(sgu.lon);
diveu.lat = nanmean(sgu.lat);
diveu.dive = d_range(1):d_range(2);
diveu.date = nanmean(sgu.date);
diveu.hour = diveu.date - fix(diveu.date);
[diveu.u10,~,~,~,~,~] = ow_2obs(diveu.lat,diveu.lon,diveu.date+10/24); % extract wind data
diveu.u10 = diveu.u10';
diveu.u10(isnan(diveu.u10)) = interp1(diveu.date(~isnan(diveu.u10)),diveu.u10(~isnan(diveu.u10)),diveu.date(isnan(diveu.u10))); % interpolate missing data
[diveu.slp,~,~,~] = nr_2obs(diveu.lat,diveu.lon,diveu.date+10/24,'pres','surface_gauss'); % Sea level pressure in Pascals
diveu.slp = diveu.slp';
diveu.units = {'lon: Degrees E';'lat: Degrees N';'dive: #';'date: Matlab date (HST)';'hour: Decimal day (HST)';'u10: wind at 10 m (m s-1)';'slp: sea level pressure (Pa)'};

% 8. Oxygen calibration on Winkler measurements
    % build comparison array
o2_comp = oxy_cal;
o2_comp.dist = NaN(height(oxy_cal),1);
o2_comp.lag = NaN(height(oxy_cal),1);
o2_comp.zdist = NaN(height(oxy_cal),1);
o2_comp.optode = NaN(height(oxy_cal),1);
o2_comp.sbe43 = NaN(height(oxy_cal),1);
for i = 1:height(oxy_cal)
    [o2_comp.lag(i),ind_t] = min(abs(diveu.date-o2_comp.date(i)));
    o2_comp.dist(i) = vdist(diveu.lat(ind_t),diveu.lon(ind_t),o2_comp.lat(i),-o2_comp.lon(i));
    [o2_comp.zdist(i),ind_z] = min(abs(sgu.depth-o2_comp.depth(i)));
    o2_comp.optode(i) = sgu.opt(ind_z,ind_t);
    o2_comp.sbe43(i) = sgu.o(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10;
ind_correction = o2_comp.dist<10000 & o2_comp.lag<0.2 & o2_comp.zdist<10 & o2_comp.depth<30;
    % Compute offset in the upper 100 m and correct oxygen from sensors
optode_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction));
sbe43_offset = nanmedian(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction));
sgu.opt = sgu.opt + optode_offset;
sgu.o = sgu.o + sbe43_offset;
    % Plot comparison results
%{
subplot(2,2,1)
scatter(o2_comp.oxygen(ind_good),o2_comp.optode(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('Optode oxygen (umol L-1)')
subplot(2,2,2)
scatter(o2_comp.oxygen(ind_good),o2_comp.sbe43(ind_good),200,o2_comp.depth(ind_good),'.')
caxis([0 500]), colormap(jet), cb = colorbar, title(cb,'depth (m)','Fontsize',16)
ylim([80 240]),xlim([80 240])
hold on, plot([80 240],[80 240],'k--'), hold off
set(gca,'Fontsize',16), xlabel('Winkler oxygen (umol L-1)'), ylabel('SBE43 oxygen (umol L-1)')
subplot(2,2,3)
hist(o2_comp.oxygen(ind_correction)-o2_comp.optode(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - optode (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
subplot(2,2,4)
hist(o2_comp.oxygen(ind_correction)-o2_comp.sbe43(ind_correction),20)
set(gca,'Fontsize',16), xlabel('Winkler - SBE43 (umol L-1)'), ylabel('n. obs.')
title('upper 100m')
%}

% 9. Fluorescence-Chlorophyll calibration from HPLC pigment measurements
% build comparison array
hplc_comp = hplc_cal(:,{'DATE','DEPTH','LAT','LON','TOTCHLA'});
hplc_comp.TPP = hplc_cal.TOTCHLA + hplc_cal.CHLC + hplc_cal.MV_CHLB + hplc_cal.FUCO + hplc_cal.HEX + hplc_cal.BUT;
hplc_comp.dist = NaN(height(hplc_cal),1);
hplc_comp.lag = NaN(height(hplc_cal),1);
hplc_comp.zdist = NaN(height(hplc_cal),1);
hplc_comp.chl1 = NaN(height(hplc_cal),1);
hplc_comp.chl2 = NaN(height(hplc_cal),1);
for i = 1:height(hplc_cal)
    [hplc_comp.lag(i),ind_t] = min(abs(diveu.date-hplc_comp.DATE(i)));
    hplc_comp.dist(i) = vdist(diveu.lat(ind_t),diveu.lon(ind_t),hplc_comp.LAT(i),-hplc_comp.LON(i));
    [hplc_comp.zdist(i),ind_z] = min(abs(sgu.depth-hplc_comp.DEPTH(i)));
    hplc_comp.chl1(i) = sgu.chl1(ind_z,ind_t);
    hplc_comp.chl2(i) = sgu.chl2(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = hplc_comp.dist<20000 & hplc_comp.lag<0.2 & hplc_comp.zdist<10;
    % Compute offset in the upper 100 m and correct oxygen from sensors
% Regression fluorescence vs. chla+chlb+chlc+psc
[param_c1,stat_c1] = robustfit(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good));
[slope_psc1,intercept_psc1,r_psc1,sslope_psc1,sintercept_psc1]  = lsqfitgm(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good));
[param_c2,stat_c2] = robustfit(hplc_comp.chl2(ind_good),hplc_comp.TOTCHLA(ind_good)+hplc_comp.TPP(ind_good));
opts = optimset('MaxIter',20000,'MaxFunEval',200000);
beta3 = lsqcurvefit('gompertzfun',[0.89 -0.31 -3000 -0.08],hplc_comp.DEPTH,hplc_comp.TOTCHLA./hplc_comp.TPP,[],[],opts);
r_sig3 = corrcoef(gompertzfun(beta3,hplc_comp.DEPTH),hplc_comp.TOTCHLA./hplc_comp.TPP);
% Calibrate fluorescence on chlorophyll a
sgu.chl1 = (param_c1(2)*sgu.chl1 + param_c1(1)).*(gompertzfun(beta3,repmat(sgu.depth,1,nd)));
sgu.chl2 = (param_c2(2)*sgu.chl2 + param_c2(1)).*(gompertzfun(beta3,repmat(sgu.depth,1,nd)));
sgu.Properties.VariableUnits({'chl1','chl2'}) = {'ng chl-a L-1','ng chl-a L-1'};
    % Plot comparison results
%{
subplot(2,2,1)
plot(hplc_comp.chl1(ind_good),hplc_comp.TPP(ind_good),'ko')
hold on, plot([nanmin(hplc_comp.chl1) nanmax(hplc_comp.chl1)],param_c1(2)*[nanmin(hplc_comp.chl1) nanmax(hplc_comp.chl1)]+param_c1(1),'k--'), hold off
xlabel('Seaglider Fluorescence 1'), ylabel('Total Photosynthetic Pigments(ng/L)'), set(gca,'Fontsize',18)
lg = legend('data','robust fit'), set(lg,'box','off','Location','NorthWest','Fontsize',18)
subplot(2,2,3)
plot(hplc_comp.chl2(ind_good),hplc_comp.TPP(ind_good),'ko')
hold on, plot([nanmin(hplc_comp.chl2) nanmax(hplc_comp.chl2)],param_c2(2)*[nanmin(hplc_comp.chl2) nanmax(hplc_comp.chl2)]+param_c2(1),'k--'), hold off
xlabel('Seaglider Fluorescence 2'), ylabel('Total Photosynthetic Pigments(ng/L)'), set(gca,'Fontsize',18)
subplot(2,2,[2 4])
plot(hplc_comp.TOTCHLA./hplc_comp.TPP,hplc_comp.DEPTH,'ko'),hold on, plot(gompertzfun(beta3,0:200),0:1:200,'k--'), hold off, set(gca,'ydir','rev','Fontsize',18)
xlabel('Chlorophyll a / Total Photosynthetic Pigments(ng/ng)'), ylabel('Depth (m)')
%}
    
% Compute characteristics on isopycnal levels (average sigma at depth bins)
    % initialize table of isopycnal data
isou = sgu;
isou.sig = sig_grid;
isou.date = NaN*isou.date; isou.hour = NaN*isou.hour; isou.lat = NaN*isou.lat; isou.lon = NaN*isou.lon;
isou.s = NaN*isou.s; isou.t = NaN*isou.t; isou.o = NaN*isou.o; isou.opt = NaN*isou.opt;
isou.chl1 = NaN*isou.chl1; isou.chl2 = NaN*isou.chl2; isou.cdom = NaN*isou.cdom;
isou.bbp470 = NaN*isou.bbp470; isou.bbp660 = NaN*isou.bbp660; isou.bbp700 = NaN*isou.bbp700;
    % interpolate at isopycnal values
for i = 1:nd
    ind_nan = isnan(sgu.sig(:,i));
    if sum(~ind_nan) ~= 0
        isou.t(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.t(~ind_nan,i),sig_grid);
        isou.s(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.s(~ind_nan,i),sig_grid);
        isou.date(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.date(~ind_nan,i),sig_grid);
        isou.hour(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.hour(~ind_nan,i),sig_grid);
        isou.o(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.o(~ind_nan,i),sig_grid);
        isou.opt(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.opt(~ind_nan,i),sig_grid);
        isou.chl1(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.chl1(~ind_nan,i),sig_grid);
        isou.chl2(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.chl2(~ind_nan,i),sig_grid);
        isou.cdom(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.cdom(~ind_nan,i),sig_grid);
        isou.bbp470(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.bbp470(~ind_nan,i),sig_grid);
        isou.bbp660(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.bbp660(~ind_nan,i),sig_grid);
        isou.bbp700(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.bbp700(~ind_nan,i),sig_grid);
    end
    clear ind_nan
end
clear sig_grid
    
% Compute surface oxygen flux
[Fd,Fp,Fc,~] = fas_L13(sgu.opt(2,:)/1000,diveu.u10,sgu.s(2,:)*35/35.16504,sgu.t(2,:),diveu.slp*0.00000986923267,'O2');
diveu.Fs = (Fd+Fp+Fc)*86400; % from s-1 to d-1
diveu.units{end+1} = 'Fs: surface oxygen flux (mol m-2 d-1)';
clear Fd Fp Fc

% InterpolateSLA on seaglider position and time
[diveu.sla,~,~,~] = av_2obs(diveu.lat,diveu.lon,diveu.date+10/24,'sla','linear');
diveu.sla = diveu.sla*100;
diveu.units{end+1} = 'sla: sea level anomaly from AVISO (cm)';

% Find mixed layer depth and density at the base of the mixed layer (0.03 kg m-3 difference from sigma at 10 m depth)
mld = sgu.sig - repmat(sgu.sig(5,:),height(sgu),1) - 0.03;
mld(mld<0) = NaN; mld(1,:) = NaN;
[sig003,ind003] = nanmin(mld);
mld003 = sgu.depth(ind003);
mld003sig = sig003 + sgu.sig(2,:) + 0.03;
diveu.mld003 = mld003';
diveu.mld003sig = mld003sig;
diveu.units{end+1} = 'mld003: mixed layer depth with drho 0.03 kg m-3 from 10m (m)';
diveu.units{end+1} = 'mld003sig: potential density anomaly at MLD (kg m-3)';
clear mld ind003 sig003 mld003 mld003sig  

%% Save new variables sgd, dived, isod,sgu, diveu and isou
    save(datafile,'sgd','dived','isod','sgu','diveu','isou')