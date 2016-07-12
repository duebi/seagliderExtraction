% extractSg148m14.m
%
% script to extract data from 2016 SeaGlider mission sg148m14
% 
% Benedetto Barone - Mar 2016

mission = 'sg148_m14';
d_range = [1 270]; 
upth = userpath; 
sgpath =  [upth(1:end-1) '/Data/seaglider/' mission];
%load([upth(1:end-1) '/Data/seaglider/ccar2015'])
%load([upth(1:end-1) '/Data/seaglider/aviso2015'])
clear upth
load([sgpath '/oxy_cal'])
load([sgpath '/chl_cal'])
%load pcpn_cal
cd(sgpath)
nd = length(d_range(1):d_range(2)); % total number of dives

%{
% 1. Collect data from all dives:
[UP,DWN,~] = sg_cat(d_range,sgpath,'userVarNames',{'vmtime','press','sigmath0','salin', ...
    'tempc','oxygen','optode_oxygen','optode_dphase_oxygen','wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'});

% 2. Rename WETLabs variables
DWN.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb650','cdom','chl2'};
UP.Properties.VariableNames({'wlbb2fl_bb1','wlbb2fl_bb2','wlbb2fl_chl','wlbbfl2_bb1','wlbbfl2_cdom','wlbbfl2_chl'}) = {'bb470','bb700','chl1','bb650','cdom','chl2'};

% 3. Clean backscattering data based on 3 std on a 20m grid (sets outliers = NaN)
for i = d_range(1):d_range(2)
    ind_d = DWN.divenum==i & ~isnan(DWN.bb470);
    [bin_d,bin,ind_good_d] = binning(DWN.bb470(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb470(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb470(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb470);
    [bin_u,bin,ind_good_u] = binning(UP.bb470(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb470(ind_u); temp_u(~ind_good_u) = NaN; UP.bb470(ind_u) = temp_u;
    clear ind_d ind_u ind_good_d ind_good_u temp_d temp_u bin bin_d bin_u
    ind_d = DWN.divenum==i & ~isnan(DWN.bb650);
    [bin_d,bin,ind_good_d] = binning(DWN.bb650(ind_d),DWN.vmdepth(ind_d),20,[10 190],'mean',3);
    temp_d = DWN.bb650(ind_d); temp_d(~ind_good_d) = NaN; DWN.bb650(ind_d) = temp_d;
    ind_u = UP.divenum==i & ~isnan(UP.bb650);
    [bin_u,bin,ind_good_u] = binning(UP.bb650(ind_u),UP.vmdepth(ind_u),20,[10 190],'mean',3);
    temp_u = UP.bb650(ind_u); temp_u(~ind_good_u) = NaN; UP.bb650(ind_u) = temp_u;
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
betasw470d = NaN(ldwn,1); betasw650d = NaN(ldwn,1); betasw700d = NaN(ldwn,1);
betasw470u = NaN(lup,1); betasw650u = NaN(lup,1); betasw700u = NaN(lup,1);
    % compute scatter at 117 deg due to seawater
for i = 1:ldwn % downcast 
    [betasw470d(i),~,~]= betasw_ZHH2009(470,DWN.tempc(i),117,DWN.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw650d(i),~,~]= betasw_ZHH2009(650,DWN.tempc(i),117,DWN.salin(i));
    [betasw700d(i),~,~]= betasw_ZHH2009(700,DWN.tempc(i),117,DWN.salin(i));
end
for i = 1:lup % upcast
    [betasw470u(i),~,~]= betasw_ZHH2009(470,UP.tempc(i),117,UP.salin(i)); % unfortunately this function doesn't work with vectors in T or S
    [betasw650u(i),~,~]= betasw_ZHH2009(650,UP.tempc(i),117,UP.salin(i));
    [betasw700u(i),~,~]= betasw_ZHH2009(700,UP.tempc(i),117,UP.salin(i));
end
    % calculation for particle backscattering coefficients
DWN.bbp470 = 2*pi*chi_p*(DWN.bb470-betasw470d);
DWN.bbp650 = 2*pi*chi_p*(DWN.bb650-betasw650d);
DWN.bbp700 = 2*pi*chi_p*(DWN.bb700-betasw700d);
UP.bbp470 = 2*pi*chi_p*(UP.bb470-betasw470u);
UP.bbp650 = 2*pi*chi_p*(UP.bb650-betasw650u);
UP.bbp700 = 2*pi*chi_p*(UP.bb700-betasw700u);
clear ldwn lup
clear chi_p betasw470d betasw470u betasw650d betasw650u betasw700d betasw700u

% 5. Save concatenated dive file
save([mission '_cat'],'UP','DWN');
%}

% 7. Put data on a regular grid
datafile = [mission(1:5) mission(7)  mission(8:9) 'data'];
depth = 2:2:1000; ld = length(depth);
[UG,DG] = sg_grid(['./' mission '_cat'],depth,'gridVar','vmdepth','diveRange',d_range(1):d_range(2),'outVars',{'lon','lat','daten','salin','tempc','sigmath0','oxygen','optode_dphase_oxygen','optode_oxygen','chl1','chl2','bbp470','bbp650','bbp700','cdom'});
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
sgd.opt = reshape(DG.optode_dphase_oxygen,ld,nd).*((sgd.sig+1000)/1000);
sgd.chl1 = reshape(DG.chl1,ld,nd);
sgd.chl2 = reshape(DG.chl2,ld,nd);
sgd.bbp470 = reshape(DG.bbp470,ld,nd);
sgd.bbp650 = reshape(DG.bbp650,ld,nd);
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
%
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

% 9. Fluorescence-Chlorophyll calibration from chl pigment measurements
% build comparison array
chl_comp = chl_cal;
chl_comp.totp = chl_cal.chla + chl_cal.phaeop;
chl_comp.dist = NaN(height(chl_cal),1);
chl_comp.lag = NaN(height(chl_cal),1);
chl_comp.zdist = NaN(height(chl_cal),1);
chl_comp.chl1 = NaN(height(chl_cal),1);
chl_comp.chl2 = NaN(height(chl_cal),1);
for i = 1:height(chl_cal)
    [chl_comp.lag(i),ind_t] = min(abs(dived.date(3:end)-chl_comp.date(i))); ind_t = ind_t +2;
    chl_comp.dist(i) = vdist(dived.lat(ind_t),dived.lon(ind_t),chl_comp.lat(i),-chl_comp.lon(i));
    [chl_comp.zdist(i),ind_z] = min(abs(sgd.depth-chl_comp.depth(i)));
    chl_comp.chl1(i) = sgd.chl1(ind_z,ind_t);
    chl_comp.chl2(i) = sgd.chl2(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = chl_comp.dist<10000 & chl_comp.lag<0.2 & chl_comp.zdist<10;
% Regression fluorescence vs. chla+phaeopigments
[param_c1,stat_c1] = robustfit(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
[slope_psc1,intercept_psc1,r_psc1,sslope_psc1,sintercept_psc1]  = lsqfitgm(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
[param_c2,stat_c2] = robustfit(chl_comp.chl2(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
opts = optimset('MaxIter',20000,'MaxFunEval',200000);
beta3 = lsqcurvefit('gompertzfun',[0.89 -0.31 -3000 -0.08],chl_comp.depth,chl_comp.chla./chl_comp.totp,[],[],opts);
r_sig3 = corrcoef(gompertzfun(beta3,chl_comp.depth),chl_comp.chla./chl_comp.totp);
% Calibrate fluorescence on chlorophyll a
sgd.chl1 = (param_c1(2)*sgd.chl1 + param_c1(1)).*(gompertzfun(beta3,repmat(sgd.depth,1,nd)));
sgd.chl2 = (param_c2(2)*sgd.chl2 + param_c2(1)).*(gompertzfun(beta3,repmat(sgd.depth,1,nd)));
sgd.Properties.VariableUnits({'chl1','chl2'}) = {'ng chl-a L-1','ng chl-a L-1'};
%{
% Plot comparison results
subplot(2,2,1)
plot(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good),'ko')
hold on, plot([nanmin(chl_comp.chl1) nanmax(chl_comp.chl1)],param_c1(2)*[nanmin(chl_comp.chl1) nanmax(chl_comp.chl1)]+param_c1(1),'k--'), hold off
xlabel('Seaglider Fluorescence 1'), ylabel('Chlorophyll-a+Phaeopigments(ng/L)'), set(gca,'Fontsize',18)
lg = legend('data','robust fit'), set(lg,'box','off','Location','NorthWest','Fontsize',18)
subplot(2,2,3)
plot(chl_comp.chl2(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good),'ko')
hold on, plot([nanmin(chl_comp.chl2) nanmax(chl_comp.chl2)],param_c2(2)*[nanmin(chl_comp.chl2) nanmax(chl_comp.chl2)]+param_c2(1),'k--'), hold off
xlabel('Seaglider Fluorescence 2'), ylabel('Chlorophyll-a+Phaeopigments(ng/L)'), set(gca,'Fontsize',18)
subplot(2,2,[2 4])
plot(chl_comp.chla./chl_comp.totp,chl_comp.depth,'ko'),hold on, plot(gompertzfun(beta3,0:200),0:1:200,'k--'), hold off, set(gca,'ydir','rev','Fontsize',18)
xlabel('Chlorophyll a / (Chlorophyll-a+Phaeopigments)(ng/ng)'), ylabel('Depth (m)')
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
isod.bbp470 = NaN*isod.bbp470; isod.bbp650 = NaN*isod.bbp650; isod.bbp700 = NaN*isod.bbp700;
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
        isod.bbp650(:,i) = interp1(sgd.sig(~ind_nan,i),sgd.bbp650(~ind_nan,i),sig_grid);
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
dived.mld003 = mld003;
dived.mld003sig = mld003sig;
dived.units{end+1} = 'mld003: mixed layer depth with drho 0.03 kg m-3 from 10m (m)';
dived.units{end+1} = 'mld003sig: potential density anomaly at MLD (kg m-3)';
clear mld ind003 sig003 mld003 mld003sig  

% UPCAST data extraction
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
sgu.opt = reshape(UG.optode_dphase_oxygen,ld,nd).*((sgu.sig+1000)/1000);
sgu.chl1 = reshape(UG.chl1,ld,nd);
sgu.chl2 = reshape(UG.chl2,ld,nd);
sgu.bbp470 = reshape(UG.bbp470,ld,nd);
sgu.bbp650 = reshape(UG.bbp650,ld,nd);
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

diveu.u10 = NaN(1,nd);
[diveu.u10(1:269),~,~,~,~,~] = ow_2obs(diveu.lat(1:269),diveu.lon(1:269),diveu.date(1:269)+10/24); % extract wind data
%diveu.u10 = diveu.u10';
diveu.u10(isnan(diveu.u10)) = interp1(diveu.date(~isnan(diveu.u10)),diveu.u10(~isnan(diveu.u10)),diveu.date(isnan(diveu.u10))); % interpolate missing data

diveu.slp = NaN(1,nd);
[diveu.slp(1:269),~,~,~] = nr_2obs(diveu.lat(1:269),diveu.lon(1:269),diveu.date(1:269)+10/24,'pres','surface_gauss'); % Sea level pressure in Pascals
%diveu.slp = diveu.slp';
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
    [o2_comp.lag(i),ind_t] = min(abs(diveu.date(3:end)-o2_comp.date(i))); ind_t = ind_t +2;
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

% 9. Fluorescence-Chlorophyll calibration from chl pigment measurements
% build comparison array
chl_comp = chl_cal;
chl_comp.totp = chl_cal.chla + chl_cal.phaeop;
chl_comp.dist = NaN(height(chl_cal),1);
chl_comp.lag = NaN(height(chl_cal),1);
chl_comp.zdist = NaN(height(chl_cal),1);
chl_comp.chl1 = NaN(height(chl_cal),1);
chl_comp.chl2 = NaN(height(chl_cal),1);
for i = 1:height(chl_cal)
    [chl_comp.lag(i),ind_t] = min(abs(diveu.date(3:end)-chl_comp.date(i))); ind_t = ind_t +2;
    chl_comp.dist(i) = vdist(diveu.lat(ind_t),diveu.lon(ind_t),chl_comp.lat(i),-chl_comp.lon(i));
    [chl_comp.zdist(i),ind_z] = min(abs(sgu.depth-chl_comp.depth(i)));
    chl_comp.chl1(i) = sgu.chl1(ind_z,ind_t);
    chl_comp.chl2(i) = sgu.chl2(ind_z,ind_t);
    clear ind_t ind_z
end
ind_good = chl_comp.dist<10000 & chl_comp.lag<0.2 & chl_comp.zdist<10;
% Regression fluorescence vs. chla+phaeopigments
[param_c1,stat_c1] = robustfit(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
[slope_psc1,intercept_psc1,r_psc1,sslope_psc1,sintercept_psc1]  = lsqfitgm(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
[param_c2,stat_c2] = robustfit(chl_comp.chl2(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good));
opts = optimset('MaxIter',20000,'MaxFunEval',200000);
beta3 = lsqcurvefit('gompertzfun',[0.89 -0.31 -3000 -0.08],chl_comp.depth,chl_comp.chla./chl_comp.totp,[],[],opts);
r_sig3 = corrcoef(gompertzfun(beta3,chl_comp.depth),chl_comp.chla./chl_comp.totp);
% Calibrate fluorescence on chlorophyll a
sgu.chl1 = (param_c1(2)*sgu.chl1 + param_c1(1)).*(gompertzfun(beta3,repmat(sgu.depth,1,nd)));
sgu.chl2 = (param_c2(2)*sgu.chl2 + param_c2(1)).*(gompertzfun(beta3,repmat(sgu.depth,1,nd)));
sgu.Properties.VariableUnits({'chl1','chl2'}) = {'ng chl-a L-1','ng chl-a L-1'};
%{
% Plot comparison results
subplot(2,2,1)
plot(chl_comp.chl1(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good),'ko')
hold on, plot([nanmin(chl_comp.chl1) nanmax(chl_comp.chl1)],param_c1(2)*[nanmin(chl_comp.chl1) nanmax(chl_comp.chl1)]+param_c1(1),'k--'), hold off
xlabel('Seaglider Fluorescence 1'), ylabel('Chlorophyll-a+Phaeopigments(ng/L)'), set(gca,'Fontsize',18)
lg = legend('data','robust fit'), set(lg,'box','off','Location','NorthWest','Fontsize',18)
subplot(2,2,3)
plot(chl_comp.chl2(ind_good),chl_comp.chla(ind_good)+chl_comp.totp(ind_good),'ko')
hold on, plot([nanmin(chl_comp.chl2) nanmax(chl_comp.chl2)],param_c2(2)*[nanmin(chl_comp.chl2) nanmax(chl_comp.chl2)]+param_c2(1),'k--'), hold off
xlabel('Seaglider Fluorescence 2'), ylabel('Chlorophyll-a+Phaeopigments(ng/L)'), set(gca,'Fontsize',18)
subplot(2,2,[2 4])
plot(chl_comp.chla./chl_comp.totp,chl_comp.depth,'ko'),hold on, plot(gompertzfun(beta3,0:200),0:1:200,'k--'), hold off, set(gca,'ydir','rev','Fontsize',18)
xlabel('Chlorophyll a / (Chlorophyll-a+Phaeopigments)(ng/ng)'), ylabel('Depth (m)')
%}

    % initialize table of isopycnal data
isou = sgu;
isou.sig = sig_grid;
isou.date = NaN*isou.date; isou.hour = NaN*isou.hour; isou.lat = NaN*isou.lat; isou.lon = NaN*isou.lon;
isou.s = NaN*isou.s; isou.t = NaN*isou.t; isou.o = NaN*isou.o; isou.opt = NaN*isou.opt;
isou.chl1 = NaN*isou.chl1; isou.chl2 = NaN*isou.chl2; isou.cdom = NaN*isou.cdom;
isou.bbp470 = NaN*isou.bbp470; isou.bbp650 = NaN*isou.bbp650; isou.bbp700 = NaN*isou.bbp700;
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
        isou.bbp650(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.bbp650(~ind_nan,i),sig_grid);
        isou.bbp700(:,i) = interp1(sgu.sig(~ind_nan,i),sgu.bbp700(~ind_nan,i),sig_grid);
    end
    clear ind_nan
end
clear sig_grid

% Compute surface oxygen flux
[Fd,Fp,Fc,~,~] = fas_L13(sgu.opt(2,:)/1000,diveu.u10,sgu.s(2,:)*35/35.16504,sgu.t(2,:),diveu.slp*0.00000986923267,'O2');
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
diveu.mld003 = mld003;
diveu.mld003sig = mld003sig;
diveu.units{end+1} = 'mld003: mixed layer depth with drho 0.03 kg m-3 from 10m (m)';
diveu.units{end+1} = 'mld003sig: potential density anomaly at MLD (kg m-3)';
clear mld ind003 sig003 mld003 mld003sig

% Save new variables sgd, dived, isod, sgu, diveu and isou
    save(datafile,'sgd','dived','isod','sgu','diveu','isou')
%%

ind_plot = 164:nd;

%figure
clf
subplot(5,2,[1 3 5])
contourf(dived.lat(ind_plot),sgd.depth,sgd.s(:,ind_plot),32.75:0.05:35.4,'edgecolor','none')
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--'), hold off
hold on, plot(dived.lat(ind_plot),dived.mld003(ind_plot),'k')
set(gca,'ydir','rev','XAxisLocation','top')
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 350])
caxis([34 35.35])
cb = colorbar; title(cb,'Salinity (g/kg)'); set(cb,'Fontsize',18)
subplot(5,2,[2 4 6])
contourf(dived.lat(ind_plot),sgd.depth,sgd.opt(:,ind_plot),0:1:235,'edgecolor','none')
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--'), hold off
hold on, plot(dived.lat(ind_plot),dived.mld003(ind_plot),'k')
set(gca,'ydir','rev','XAxisLocation','top')
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 350])
caxis([140 226])
cb = colorbar; title(cb,'Oxygen (\mumol/L)'); set(cb,'Fontsize',18)
subplot(5,2,[7 9])
contourf(dived.lat(ind_plot),sgd.depth,sgd.chl2(:,ind_plot),0:0.01:3,'edgecolor','none','Linewidth',1.5)
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--','Linewidth',1.5), hold off
hold on, plot(dived.lat(ind_plot),dived.mld003(ind_plot),'k')
set(gca,'ydir','rev','Fontsize',18)
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 200])
caxis([0 1.1])
cb = colorbar; title(cb,'Chlorophyll-a (mg m-3)'); set(cb,'Fontsize',18)
subplot(5,2,[8 10])
contourf(dived.lat(ind_plot),sgd.depth,sgd.bbp650(:,ind_plot),0:2.5e-5:2e-3,'edgecolor','none')
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--','Linewidth',1.5), hold off
hold on, plot(dived.lat(ind_plot),dived.mld003(ind_plot),'k','Linewidth',1.5)
set(gca,'ydir','rev','Fontsize',18)
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 200])
caxis([0.05e-3 0.55e-3])
cb = colorbar; title(cb,'bbp650 (1/m)'); set(cb,'Fontsize',18)
colormap(jet(50))

%%

ind_plot = 164:nd;

%figure
clf
subplot(3,3,[1 4 7])
% coastline
load hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii

% Aviso altimetry for April 30
aviso_filename = ['/Users/benedetto/Documents/MATLAB/Data/km1605/altimetry/aviso/nrt_global_allsat_msla_h_20160430_20160506.nc'];
sla.lat = ncread(aviso_filename,'lat');
sla.lon = ncread(aviso_filename,'lon');
i_lat = find(sla.lat>=20 & sla.lat<=25.5);
nlat = length(i_lat);
i_lon = find(sla.lon>=200.5 & sla.lon<=206);
nlon = length(i_lon);
sla.lat = sla.lat(i_lat);
sla.lon = -(sla.lon(i_lon)-360);
sla.sla = 100*ncread(aviso_filename,'sla',[i_lon(1) i_lat(1) 1],[nlon nlat 1])';

hold on, plot(-dived.lon(ind_plot),dived.lat(ind_plot),'k-','Linewidth',2), hold off

hold on
    % Coastline
    for i = 1:(length(ind_n)-1)
        patch(coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none','FaceAlpha',1)
    end
    % Station ALOHA
     viscircles([158 22.75],6/59.79,'edgecolor','k','linewidth',2,'DrawBackgroundCircle',0);
     text(158.2,22.95,'  Sta. ALOHA','color','k','Fontsize',18);
     % Altimetry
     [c,h] = contour(sla.lon,sla.lat,sla.sla,'k--');
     clabel(c,h,'Fontsize',16);
hold off
axis equal
xlabel('Longitude W')
ylabel('Latitude N')
set(gca,'xdir','rev','Fontsize',18)
xlim([156 158.5]), ylim([21 25])
title('Sea level anomaly (cm) & Seaglider track','FontWeight','normal')

subplot(3,3,[2 3])
contourf(dived.lat(ind_plot),sgd.depth,sgd.s(:,ind_plot),32.75:0.05:35.4,'edgecolor','none')
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--'), hold off
hold on, plot(dived.lat(ind_plot),mld003(ind_plot),'k')
set(gca,'ydir','rev','XAxisLocation','top')
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 200])
caxis([34 35.35])
cb = colorbar; title(cb,'Salinity (g/kg)'); set(cb,'Fontsize',18)
subplot(3,3,[5 6])
contourf(dived.lat(ind_plot),sgd.depth,sgd.chl2(:,ind_plot),0:0.01:3,'edgecolor','none','Linewidth',1.5)
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--','Linewidth',1.5), hold off
hold on, plot(dived.lat(ind_plot),mld003(ind_plot),'k'), hold off
%hold on, plot(dived.lat(ind_plot),mld003(ind_plot),'k')
set(gca,'ydir','rev','Fontsize',18)
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 200])
caxis([0 1.1])
cb = colorbar; title(cb,'Chlorophyll-a (mg m-3)'); set(cb,'Fontsize',18)
subplot(3,3,[8 9])
contourf(dived.lat(ind_plot),sgd.depth,sgd.bbp650(:,ind_plot),0:2.5e-5:2e-3,'edgecolor','none')
hold on, contour(dived.lat(ind_plot),sgd.depth,sgd.sig(:,ind_plot),[23.5 24.75 25.5],'k--','Linewidth',1.5), hold off
hold on, plot(dived.lat(ind_plot),mld003(ind_plot),'k'), hold off
%hold on, plot(dived.lat(ind_plot),mld003(ind_plot),'k','Linewidth',1.5)
set(gca,'ydir','rev','Fontsize',18)
ylabel('Depth (m)')
xlabel('Latitude N')
%datetick('x','mm/dd'), xlabel('Month/day 2016 (HST)')
xlim([dived.lat(ind_plot(1)) dived.lat(ind_plot(end))])
ylim([0 200])
caxis([0.05e-3 0.55e-3])
cb = colorbar; title(cb,'bbp650 (1/m)'); set(cb,'Fontsize',18)
colormap(jet(50))

%%
ind_dives = 164:nd;
%ind_dives = 89:123;

figure
clf
% coastline
load hawaii.dat
coastline = hawaii;
coastline(:,1) = -coastline(:,1);
ind_n = find(isnan(coastline(:,1)));
clear hawaii

hold on, plot(-dived.lon(ind_dives),dived.lat(ind_dives),'k--','Linewidth',3), hold off
%{
hold on, surf(-sgd.lon,sgd.lat,repmat(sgd.depth,1,nd),sgd.chl2,'edgecolor','none'),hold off
set(gca,'xdir','rev','zdir','rev','zlim',[0 200],'Color','none','Fontsize',18)
hold on,plot3(-dived.lon,dived.lat,mld003,'w','Linewidth',3), hold off
caxis([0 7.5])
colormap(jet)
shading interp
%}
hold on, surf(-sgd.lon(:,ind_dives),sgd.lat(:,ind_dives),repmat(sgd.depth,1,length(ind_dives)),sgd.s(:,ind_dives),'edgecolor','none'),hold off
set(gca,'xdir','rev','zdir','rev','zlim',[0 200],'Color','none','Fontsize',18)
hold on,plot3(-dived.lon(ind_dives),dived.lat(ind_dives),mld003(ind_dives),'w','Linewidth',3), hold off
caxis([34.35 35.4])
colormap(jet)
shading interp


hold on
    % Coastline
    for i = 1:(length(ind_n)-1)
        patch(coastline((ind_n(i)+1):(ind_n(i+1)-1),1),coastline((ind_n(i)+1):(ind_n(i+1)-1),2),'k','edgecolor','none','FaceAlpha',1)
    end
    % Station ALOHA
    plot(158,22.75,'k.','Markersize',30)
    %hcirc = viscircles([158 22.75],6/59.79,'edgecolor','k','linewidth',2);
    text(158,22.75,'  Sta. ALOHA','color','k','Fontsize',16);
hold off

set(gca,'DataAspectRatio',[1 1 120],'xlim',[156 158.5],'ylim',[21.25 24.5])
grid on

latlim = [20.5 24.5];
lonlim = [201 205];
% Download aviso altimetry
dtemp = datevec(datenum(clock) +10/24);
today_str = num2str(dtemp(1)*10000 + dtemp(2)*100  + dtemp(3));
aviso_filename = ['../../km1605/altimetry/aviso/nrt_global_allsat_msla_h_' today_str '_' today_str '.nc'];
sla.lat = ncread(aviso_filename,'lat');
sla.lon = ncread(aviso_filename,'lon');
i_lat = find(sla.lat>=latlim(1)-0.5 & sla.lat<=latlim(2)+0.5);
nlat = length(i_lat);
i_lon = find(sla.lon>=lonlim(1)-0.5 & sla.lon<=lonlim(2)+0.5);
nlon = length(i_lon);
sla.lat = sla.lat(i_lat);
sla.lon = -(sla.lon(i_lon)-360);
sla.sla = 100*ncread(aviso_filename,'sla',[i_lon(1) i_lat(1) 1],[nlon nlat 1])';

[~,c_lon] = max(max(sla.sla));
[~,c_lat] = max(max(sla.sla'));
maxsla.lon = sla.lon(c_lon);maxsla.lat = sla.lat(c_lat);
hold on, [c,h] = contour(sla.lon,sla.lat,sla.sla,-32:4:32,'edgecolor','k'); hold off
clabel(c,h,'Fontsize',18)

xlabel('Longitude (deg W)'), ylabel('Latitude (deg N)'),zlabel('Depth (m)')