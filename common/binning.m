function [bin_data,bin,ind_good] = binning(rawdata,pressure,binwidth,binrange,method,nstd)
%
% Simple binning of data with or without outlier removal. Data are binned
% independently columnwise.
%
% [bin_data,bin,ind_good] = binning(rawdata,pressure,binwidth,binrange,method,nstd)
%
% INPUT:    rawdata: matrix of raw data (n*m)
%           pressure: raw pressure values (n*1)
%           binwidth: size of the pressure bin (DEFAULT: 2)
%           binrange: first and last bin center (DEFAULT: [min(pressure)+binwidth/2 max(pressure)])
%           method: can be 'mean' or 'median'
%           nstd =  number of standard deviation from the mean to define
%           outliers (DEFAULT: 'no' outliers removal)
% OUTPUT:   bin_data: data after binning (nbins*m)
%           bin: central pressure value at bins (nbins*1)
%           ind_good: logical indeces of good values(n*m)
%
% Autohr: Benedetto Barone

% Set default values for bin range and binwidth
if nargin < 3
    binwidth = 2;
end
if nargin < 4
    binrange = [min(pressure)+binwidth/2 max(pressure)];
end
if nargin < 5
    method = 'mean';
end
if nargin < 6
    nstd = 'no';
end

%Check size
lp = length(pressure);
size_raw = size(rawdata);
m = size_raw(2);
n = size_raw(1);
if n ~= lp
    error('Incorrect number of rows','the number of rows of rawdata and pressure must be the same')
end
% Bin parameters
bin = binrange(1):binwidth:binrange(2); bin = bin';
lb = length(bin);
% Initialize bin variables
bin_data = NaN(lb,m);
ind_good = logical(ones(lb,m));

% ---------
% Bin data 
% ---------
% no outliers removal
if strcmp(nstd,'no') == 1
    for i = 1:lb % Go through all the bins
        % Isolate data from the bin
        ind_bin = pressure < bin(i)+binwidth/2 & pressure >= bin(i)-binwidth/2;
        data_bin = rawdata(ind_bin,:);
        bin_data(i,:) = eval(['nan' method '(data_bin)']);
    end
    % with outliers removal
else
    for i = 1:lb % Go through all the bins
        % Isolate data from the bin
        ind_bin = pressure < bin(i)+binwidth/2 & pressure >= bin(i)-binwidth/2;
        data_bin = rawdata(ind_bin,:);
        for j = 1:m
            mean_bin = nanmean(data_bin(:,j));
            std_bin = nanstd(data_bin(:,j));
            good_bin = data_bin(:,j) <= mean_bin + nstd*std_bin & data_bin(:,j) >= mean_bin - nstd*std_bin;
            ind_good(ind_bin,j) = good_bin;
            bin_data(i,j) = eval(['nan' method '(data_bin(good_bin,j))']);
            clear mean_bin std_bin good_bin
        end
    end
end