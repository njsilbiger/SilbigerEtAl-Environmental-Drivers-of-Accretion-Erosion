%% %Hawaiian Archipelago data analysis
%created by Nyssa Silbiger 1/6/15
%updated 3/16/17


%calculate relationships between accretion and erosion and different
%environmental variables (chemistry,Fish, benthic) along the HA

%% %Load accretion-erosion data
%accretion erosion and meta data
load('HAAccretionErosion.mat');

%% %Load chemistry data

%[~, ~, raw] = xlsread('\Chemistry_MATLAB.xlsx','Sheet1','C2:W59');
[~, ~, raw] = xlsread('\Chemistry_MATLAB.xlsx','Revisions','C2:W63');
%[~, ~, raw] = xlsread('\Chemistry_MATLAB.xlsx','Original','C2:W66');

raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[2,3,4]);
raw = raw(:,[1,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));

% Create dataset array
Chem = dataset;

% Allocate imported array to column variable names
Chem.Year = data(:,1);
Chem.LocationCo = cellVectors(:,1);
Chem.Location = cellVectors(:,2);
Chem.REASite = cellVectors(:,3);
Chem.Latitude = data(:,2);
Chem.Longitude = data(:,3);
Chem.UTCDateTim = data(:,4);
Chem.SampleDep = data(:,5);
Chem.PChl = data(:,6);
Chem.PO = data(:,7);
Chem.Si = data(:,8);
Chem.NO30 = data(:,9);
Chem.NO2 = data(:,10);
Chem.NO32 = data(:,11); %nitrate+nitrite
Chem.DIC = data(:,12);
Chem.TA = data(:,13);
Chem.Salinity4 = data(:,14);
Chem.TAnorm = data(:,15); %salinity normalized TA
Chem.Temp = data(:,16);
Chem.Press = data(:,17);
Chem.DICnorm=data(:,18); %salinity noramalized DIC
% Clear temporary variables
clearvars data raw cellVectors R;

%% %Calculate Carb Chem Params

[CO2Data,CO2_Headers,NICEHEADERS]=CO2SYS(Chem.TA,Chem.DIC,1,2,Chem.Salinity4,Chem.Temp,Chem.Temp,0,0,0,0,1,10,1);

%Add the CC params to the Chem dataset
Chem.Omega=CO2Data(:,31); 
Chem.pH=CO2Data(:,37);
Chem.HCO3=CO2Data(:,21); 
Chem.CO2=CO2Data(:,23); 
Chem.CO3=CO2Data(:,22); 

%Only use benthic samples
a=find(Chem.SampleDep<1.5); %only use date >1.5 meters deep
Chem(a,:)=[];


%% Load Fish Data
filename = '\Fish_mid_buffered.csv';
delimiter = ',';
startRow = 2;

% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

% Open the text file.
fileID = fopen(filename,'r');

dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[6,7,8,9,10,11,12,14,15,17,21,22,23,24,25,27,28,29]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [6,7,8,9,10,11,12,14,15,17,21,22,23,24,25,27,28,29]);
rawCellColumns = raw(:, [1,2,3,4,5,13,16,18,19,20,26]);


% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

% Create output variable
Fish = dataset;
Fish.TAXON = rawCellColumns(:, 1);
Fish.METHOD = rawCellColumns(:, 2);
Fish.SITE = rawCellColumns(:, 3);
Fish.DATE = rawCellColumns(:, 4);
Fish.REEF_ZONE = rawCellColumns(:, 5);
Fish.LATITUDE = cell2mat(rawNumericColumns(:, 1));
Fish.LONGITUDE = cell2mat(rawNumericColumns(:, 2));
Fish.HARD_CORAL = cell2mat(rawNumericColumns(:, 3));
Fish.MA = cell2mat(rawNumericColumns(:, 4));
Fish.TA = cell2mat(rawNumericColumns(:, 5));
Fish.CCA = cell2mat(rawNumericColumns(:, 6));
Fish.SAND = cell2mat(rawNumericColumns(:, 7));
Fish.SPECIES = rawCellColumns(:, 6);
Fish.ABUND_M2 = cell2mat(rawNumericColumns(:, 8));
Fish.BIO_M2 = cell2mat(rawNumericColumns(:, 9));
Fish.FAMILY = rawCellColumns(:, 7);
Fish.X = cell2mat(rawNumericColumns(:, 10));
Fish.Trophic5 = rawCellColumns(:, 8);
Fish.FAMILY_1 = rawCellColumns(:, 9);
Fish.variable = rawCellColumns(:, 10);
Fish.Biomass_corr = cell2mat(rawNumericColumns(:, 11));
Fish.Abund_corr = cell2mat(rawNumericColumns(:, 12));
Fish.meanSHcalc = cell2mat(rawNumericColumns(:, 13));
Fish.DEPTH = cell2mat(rawNumericColumns(:, 14));
Fish.YEAR = cell2mat(rawNumericColumns(:, 15));
Fish.REASite = rawCellColumns(:, 11);
Fish.Latitude_1 = cell2mat(rawNumericColumns(:, 16));
Fish.Longitud_1 = cell2mat(rawNumericColumns(:, 17));
Fish.Distance = cell2mat(rawNumericColumns(:, 18));
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

%%  Import SST Satallite Data
[~, ~, raw] = xlsread('\SSTDataAll.xlsx','Sheet1','A2:E30');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,1);
raw = raw(:,[2,3,4,5]);

% Create output variable
data = reshape([raw{:}],size(raw));

% Create table
SST = table;

% Allocate imported array to column variable names
SST.Site = cellVectors(:,1);
SST.SSTMean = data(:,1);
SST.SSTMax = data(:,2);
SST.SSTMin = data(:,3);
SST.SSTVar = sqrt(data(:,4)); %make std instead of var

% Clear temporary variables
clearvars data raw cellVectors;
%% %take averages of all chem data
[Site,PO,PO_std]=grpstats(Chem.PO,Chem.REASite,{'gname','mean','std'});
[Si,Si_std]=grpstats(Chem.Si,Chem.REASite,{'mean','std'});
[NO32,NO32_std]=grpstats(Chem.NO32,Chem.REASite,{'mean','std'});
[DIC,DIC_std]=grpstats(Chem.DICnorm,Chem.REASite,{'mean','std'}); %DIC normalized
[TA,TA_std]=grpstats(Chem.TAnorm,Chem.REASite,{'mean','std'}); %TA normalized
Lat=grpstats(Chem.Latitude,Chem.REASite,'mean');
[Omega,Omega_std]=grpstats(Chem.Omega,Chem.REASite,{'mean','std'});
[pH,pH_std]=grpstats(Chem.pH,Chem.REASite,{'mean','std'});
[HCO3,HCO3_std]=grpstats(Chem.HCO3,Chem.REASite,{'mean','std'});
[CO3,CO3_std]=grpstats(Chem.CO3,Chem.REASite,{'mean','std'});
[CO2,CO2_std]=grpstats(Chem.CO2,Chem.REASite,{'mean','std'});
Depth=grpstats(Chem.SampleDep,Chem.REASite,'mean');
NO2=grpstats(Chem.NO2,Chem.REASite,'mean');

%create chem dataset with averages at site level
ChemNew=dataset;
ChemNew.Si=Si; ChemNew.Si_std=Si_std;
ChemNew.NO32=NO32; ChemNew.NO32_std=NO32_std;
ChemNew.NO2=NO2;
ChemNew.PO=PO; ChemNew.PO_std=PO_std;
ChemNew.DIC=DIC; ChemNew.DIC_std=DIC_std;
ChemNew.TA=TA; ChemNew.TA_std=TA_std;
ChemNew.Omega=Omega; ChemNew.Omega_std=Omega_std;
ChemNew.pH=pH; ChemNew.pH_std=pH_std;
ChemNew.HCO3=HCO3; ChemNew.HCO3_std=HCO3_std;
ChemNew.CO2=CO2; ChemNew.CO2_std=CO2_std;
ChemNew.CO3=CO3; ChemNew.CO3_std=CO3_std;
ChemNew.Site=Site;
ChemNew.Depth=Depth;
ChemNew.Lat=Lat;
%ChemNew.Site=Site;
%export(ChemNew,'XLSfile','AllChem.xlsx')
%export a table with all the chem data 
%% %Accretion and Erosion Means
%normalize Percent net change so that it is percent change per year

%Remove the bad ct scan data (these all had problems with the scans and are unusable)
no=[636,642,1004,697,685,686,677,665,661,633,688];
nums=zeros(11,1);
for i=1:length(nums)
    nums(i)=find(Data.BlockID==no(i));
end

Data(nums,:)=[];

Data.PercentNetChange=Data.PercentNetChange.*(Data.Days./365);
[AccMean_S,AccSEM_S,n_S,name_S]=grpstats(Data.AccRate, Data.Site, {'mean',  'sem','numel','gname',});

[ErMean_S,ErSEM_S]=grpstats(Data.ErRate, Data.Site, {'mean', 'sem'});
[TotalMean_S,TotalSEM_S]=grpstats(Data.PercentNetChange, Data.Site, {'mean', 'sem'});
%create accretion-erosion dataset at site level
AccErMean=dataset;
AccErMean.Site=name_S; 
AccErMean.Accretion_mean=AccMean_S;AccErMean.Er_sem=ErSEM_S;
AccErMean.Erosion_mean=ErMean_S;AccErMean.Acc_sem=AccSEM_S;
AccErMean.TotalMean=TotalMean_S;
AccErMean.TotalSEM=TotalSEM_S;

%% %Analyze fish data
%Fish data from 2008 - 2014
%Weight the biomasses by distance from site

%Linear Weight
Weight=1./(ones(length(Fish.Distance),1)+Fish.Distance); %higher distances have less weight
Fish.WeightedBiomass=Weight.*Fish.Biomass_corr; %multiply biomass by weight

%Calcualte average fish biomass for all herbivores
c=find(strcmp(Fish.Trophic5,'H')); %find sherbovires
Fish_Herb=Fish(c,:); %re-create dataset with only Herbovores

[Site,Herb_Mean,Herb_sem,n]=grpstats(Fish_Herb.WeightedBiomass,Fish_Herb.REASite,{'gname','mean','sem','numel'});

%show all the data by site and year
[IDh,HerbM,HerbN,Herbsem]=grpstats(Fish_Herb.WeightedBiomass,{Fish_Herb.REASite,Fish_Herb.YEAR},{'gname','mean','numel','sem'});

%create dataset for herbivores
Herbs=dataset;
Herbs.Site=Site;
Herbs.Herb_Mean=Herb_Mean;
Herbs.Herb_sem=Herb_sem;
%Herbs.N=n;
b=find(strcmp(Herbs.Site,'MID-H11')); %remove MID H11
Herbs(b,:)=[];

%% %Load LPI data

[~, ~, raw, dateNums] = xlsread('BenthicAll.xlsx','BenthicAll','A2:U1035','',@convertSpreadsheetDates);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[3,5,6,10,11,12,14,18]);
raw = raw(:,[1,2,4,7,8,9,13,15,16,17,19,20,21]);
dateNums = dateNums(:,[1,2,4,7,8,9,13,15,16,17,19,20,21]);

% Replace date strings by MATLAB serial date numbers (datenum)
R = ~cellfun(@isequalwithequalnans,dateNums,raw) & cellfun('isclass',raw,'char'); % Find spreadsheet dates
raw(R) = dateNums(R);

% Create output variable
data = reshape([raw{:}],size(raw));

% Create dataset array
LPI = dataset;

% Allocate imported array to column variable names
LPI.OBJECTID = data(:,1);
LPI.SurveyDate = data(:,2);
LPI.Site = cellVectors(:,1);
LPI.Transect = data(:,3);
LPI.Region = cellVectors(:,2);
LPI.Island = cellVectors(:,3);
LPI.Latitude = data(:,4);
LPI.Longitude = data(:,5);
LPI.Depth = data(:,6);
LPI.BenthicCode = cellVectors(:,4);
LPI.Taxon = cellVectors(:,5);
LPI.SampMethod = cellVectors(:,6);
LPI.PctCover = data(:,7);
LPI.ReefZone = cellVectors(:,7);
LPI.UniqueSurv = data(:,8);
LPI.UniqueRepI = data(:,9);
LPI.OBJECTID_1 = data(:,10);
LPI.REASite = cellVectors(:,8);
LPI.Latitude1 = data(:,11);
LPI.Longitude1 = data(:,12);
LPI.Year = data(:,13);

% Clear temporary variables
clearvars data raw dateNums cellVectors R;

%% Species codes
% Import the data
[~, ~, raw] = xlsread('LPI_Mary.xlsx','SpeciesTable','A2:F221');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[1,2,3,4,5,6]);

% Create dataset array
Species = dataset;

% Allocate imported array to column variable names
Species.BenthicCode = cellVectors(:,1);
Species.KINGDOM = cellVectors(:,2);
Species.Taxon = cellVectors(:,3);
Species.BENTHICCATEGORY = cellVectors(:,4);
Species.COMMON_NAME = cellVectors(:,5);
Species.CATEGORY_GENUS = cellVectors(:,6);

% Clear temporary variables
clearvars data raw cellVectors;

%% join data by code 
%this gives all the data a benthic category so that I can calculate percent
%cover
%cover
%data=join(LPI,Species,'Key','BenthicCode','Type','outer');
data=join(LPI,Species,'Key','BenthicCode','Type','outer');

data.Properties.VarNames{10}='BenthicCode';data.Properties.VarNames{11}='Taxon';
ix=ismissing(data);%remove rows with missing data
data= data(~any(ix,2),:);
%remove unnecessary vars
data.BenthicCode_right=[]; data.Taxon_right=[];

%calculate percent cover by benthic category
[ID,cover]=grpstats(data.PctCover,{data.UniqueRepI,data.REASite,data.Year,data.BENTHICCATEGORY},{'gname','sum'});

[ID2,num]=grpstats(data.UniqueRepI,{data.REASite,data.Year},{'gname','numel'});


%create a dataset
PercentCover=dataset;
PercentCover.Survey=ID(:,1);
PercentCover.Site=ID(:,2);
PercentCover.Year=ID(:,3);
PercentCover.Code=ID(:,4);
PercentCover.Cover=cover;

%add replace cyan, dead, and inverts with other so that I can get the
%correct means and SEM

PercentCover.Code = strrep(PercentCover.Code, 'CYAN', 'OTHER');
PercentCover.Code = strrep(PercentCover.Code, 'DEAD', 'OTHER');
PercentCover.Code = strrep(PercentCover.Code, 'INVT', 'OTHER');

%show all the data by site and year
[IDm,coverm,coverN,coversem]=grpstats(PercentCover.Cover,{PercentCover.Site,PercentCover.Year},{'gname','mean','numel','sem'});

%organize it like a pivot table with average percent cover across year
PC=dataset2cell(PercentCover);
PC(1,:)=[];
BenthicCover=pivottable(PC, 2,4,5, @mean);
VarNames=BenthicCover(1,:);
BenthicCover(1,:)=[];
%replace empty cells with zeros
emptyIndex = cellfun(@isempty,BenthicCover);       %# Find indices of empty cells
BenthicCover(emptyIndex) = {0};                    %# Fill empty cells with 0

%convert back to dataset
BenthicCover=[VarNames;BenthicCover];
BenthicCover=cell2dataset(BenthicCover);
BenthicCover.Properties.VarNames{1}='Site';

%remove the uknown column... there is 0%
BenthicCover.UKN=[];

%add cynan, dead and invt together to get a proper sem
%BenthicCover.OTHER=BenthicCover.CYAN+BenthicCover.DEAD+BenthicCover.INVT;

%calculate the st
BenthicCoversd=pivottable(PC, 2,4,5, @std);
VarNames=BenthicCoversd(1,:);
BenthicCoversd(1,:)=[];
%replace empty cells with zeros
emptyIndex = cellfun(@isempty,BenthicCoversd);       %# Find indices of empty cells
BenthicCoversd(emptyIndex) = {0};                    %# Fill empty cells with 0
%convert back to dataset
BenthicCoversd=[VarNames;BenthicCoversd];
BenthicCoversd=cell2dataset(BenthicCoversd);
BenthicCovesdr.Properties.VarNames{1}='Site';
%remove the uknown column... there is 0%
BenthicCoversd.UKN=[];


%calculate how many samples used to make the mean
BenthicCoverN=pivottable(PC, 2,4,5, @(x)(numel(x)));
VarNames=BenthicCoverN(1,:);
BenthicCoverN(1,:)=[];
%replace empty cells with zeros
emptyIndex = cellfun(@isempty,BenthicCoverN);       %# Find indices of empty cells
BenthicCoverN(emptyIndex) = {0};                    %# Fill empty cells with 0
%convert back to dataset
BenthicCoverN=[VarNames;BenthicCoverN];
BenthicCoverN=cell2dataset(BenthicCoverN);
BenthicCoverN.Properties.VarNames{1}='Site';
%remove the uknown column... there is 0%
BenthicCoverN.UKN=[];

%find max number or reps in each row because this is really the number of
%transects... this is because some transects didn't have some species so
%instead of CRED putting a zero they did list the species.  I need this to
%calculate the correct SEM for benthic data
maxn=repmat(max(double(BenthicCoverN(:,2:end)),[],2),[1,6]);

%here is the SEM 
BenthicSEM=double(BenthicCoversd(:,2:end))./sqrt(maxn);

% I also need to add the 3 columns that are %other here because I need to
% correct means and sd. 

%  communities by genus
[ID2,cover2]=grpstats(data.PctCover,{data.UniqueRepI,data.REASite,data.Year,data.CATEGORY_GENUS},{'gname','sum'});

%create a dataset
PercentCover2=dataset;
PercentCover2.Survey=ID2(:,1);
PercentCover2.Site=ID2(:,2);
PercentCover2.Year=ID2(:,3);
PercentCover2.Code=ID2(:,4);
PercentCover2.Cover=cover2;

%organize it like a pivot table with average percent cover across year
PC2=dataset2cell(PercentCover2);
PC2(1,:)=[];
BenthicCover2=pivottable(PC2, 2,4,5, @mean);
VarNames=BenthicCover2(1,:);
BenthicCover2(1,:)=[];
%replace empty cells with zeros
emptyIndex = cellfun(@isempty,BenthicCover2);       %# Find indices of empty cells
BenthicCover2(emptyIndex) = {0};                    %# Fill empty cells with 0

%convert back to dataset
BenthicCover2=[VarNames;BenthicCover2];
BenthicCover2=cell2dataset(BenthicCover2);
BenthicCover2.Properties.VarNames{1}='Site';

%remove the uknown column... there is 0%
BenthicCover2.UKN=[];

%clear unnecessary vars
clear VarNames emptyIndex Spcies PC LPI ID cover

%% %Join all data

AllData=join(AccErMean,ChemNew); % accretion erosion and chemistry
AllData=join(AllData,Herbs); % add the herbivores
AllData=join(AllData,BenthicCover); %Benthic by code
AllData=dataset2table(AllData); %conver to table array for R2014b version
AllData=join(AllData,SST); %join the SST data


%Data2=Data;

%% join all data into one dataset

AllData2=AllData;

% %Bring in the Wave energy data
load('wave2')
wave=dataset2table(wave);
AllData2=join(AllData2,wave);

%% clear all unnecessary variables
clearvars -except name_S n_S wave AllData2 AllData AccerMean AccerMean2 Chem ChemNew Data Data2 Fish Fish_Herb Fish_Scarides Herbs Scarides BenthicCover

%sort data by site
 [i,j]=sort(AllData2.Site);
 AllData2=AllData2(j,:);
 
 FFS=1:5;
 KUR=6:10;
 LIS=11:15;
 MAI=16:19;
 OAH=20:24;
 PHR=24:29;
 
  %% Normalize Data and group bio bio, chem, and physical parameters
  %first remove a bunch of the unnecessary variables that I dont want in the
 %analysis
remove=[2:7,9,11,12,14:2:30];

% extract the accretion and erosion data
 AccEr=AllData2(:,2:7);

 AllData3=AllData2;
 AllData3(:,remove)=[]; 

%replace Na with place holders... averages of the the columns (these are
%the 3 daya points missing from MauiA27).  Can't run the model selection
%with an NA.  
AllData3.Si(19)=1.3433; AllData3.NO32(19)=0.346; AllData3.PO(19)=0.063;

 %% Biology
%BioPCA=AllData3(:,13:21);
BioPCA=table(log(AllData3.Herb_Mean), log(AllData3.CORL+1),log(AllData3.CALG+1),...
    log(AllData3.MALG+1),log(AllData3.TALG+1), log(AllData3.SAND+1),...
    log(AllData3.OTHER+1));
BioPCA.Properties.VariableNames={'Herb_Mean','Coral','CAlg','MAlg','TAlg','Sand','Other'};

%Physics
PhysPCA=table(AllData3.Depth, AllData3.SSTMean, AllData3.SSTMax, AllData3.SSTMin,...
    log(AllData3.SSTVar), log(AllData3.Energy_mean), log(AllData3.Energy_max),...
    AllData3.Energy_sum, log(AllData3.Energy_std));

PhysPCA.Properties.VariableNames={'Depth','SSTMean','SSTMax','SSTMin','SSTStd',...
    'Energy_mean','Energy_max','Energy_sum','Energy_std'};

%Chemistry
ChemPCA=table(log(AllData3.Si),log(AllData3.NO32),log(AllData3.PO),...
    log(AllData3.Omega), log(AllData3.pH), log(AllData3.TA), log(AllData3.DIC));
ChemPCA.Properties.VariableNames={'Si','NO32','PO','Omega','pH','TA','DIC'};

%% run analysis with all data, not just site means
 %repeat environmental data so that I can run models without using site
 %averages for accretion-erosion data
  [k,ix]=sort(name_S);
 name_S2=name_S(ix); n_S2=n_S(ix); %sort the data by site name

 %repeat the enviro data correct number of times 
 index = zeros(1, sum(n_S2));
index([1; cumsum(n_S2(1:end-1,:))+1]) = 1;
 AllData4=AllData3(cumsum(index), :);
 
 %sort accretion erosion data so that it matches the enviro data
 [k,ix]=sort(Data.Site);
 Data2=Data(ix,:);
 
 %% Put all the normalized chemistry data together
 All=[ChemPCA,BioPCA,PhysPCA];
 %% %Run NWHI and MHI separately
NWHI=[1:15,25:29]; %The NWHI sites for averages
MHI=[16:24]; %the MHI sites for averages

NWHIa=[1:60,104:122]; %NWHI sites for all data
MHIa=[61:103]; % MHI sites for all data

%% Run the model selection
%ranking all the AIC values:
%HA All for accretion
X=All(cumsum(index),:); %all environmental variables
y=log(Data2.AccRate);

names={'Si','N+N','PO_4^{3-}','\Omega_{arag}','pH','TA','DIC',...
    'Herb','%Coral','%Calc','%Macroalg','%Turf','%Sand',...
    '%Other','Depth','mean SST','max SST','min SST','std SST',...
    'mean Energy','max Energy ','sum Energy ','std Energy'};


for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

AccretionAll=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s', p');
AccretionAll.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope', 'p'};
[d,si]=sort(-AccretionAll.DeltaAIC);
AccretionAll=AccretionAll(si,:);
AccretionAll.Rank=si; % this put the data in same order every time for easy comparision
AccretionAll.L=exp(0.5*AccretionAll.DeltaAIC);
AccretionAll.Weight=AccretionAll.L./sum(AccretionAll.L);
%HA all for erosion----------------------------
y=log(Data2.ErRate);

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

ErosionAll=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s', p');
ErosionAll.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope' ,'p'};
[d,si]=sort(-ErosionAll.DeltaAIC);
ErosionAll=ErosionAll(si,:);
ErosionAll.Rank=si;
ErosionAll.L=exp(0.5*ErosionAll.DeltaAIC);
ErosionAll.Weight=ErosionAll.L./sum(ErosionAll.L);

%HA all for net---------------------------------

y=(Data2.PercentNetChange);

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

NetAll=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s', p');
NetAll.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-NetAll.DeltaAIC);
NetAll=NetAll(si,:);
NetAll.Rank=si;
NetAll.L=exp(0.5*NetAll.DeltaAIC);
NetAll.Weight=NetAll.L./sum(NetAll.L);

% MHI all Accretion------------------------------
X=X(MHIa,:); %all environmental variables
y=log(Data2.AccRate(MHIa));

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

AccretionAllM=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s',p');
AccretionAllM.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-AccretionAllM.DeltaAIC);
AccretionAllM=AccretionAllM(si,:);
AccretionAllM.Rank=si;
AccretionAllM.L=exp(0.5*AccretionAllM.DeltaAIC);
AccretionAllM.Weight=AccretionAllM.L./sum(AccretionAllM.L);

%MHI Erosion----------------------------------------------
y=log(Data2.ErRate(MHIa));

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

ErosionAllM=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s',p');
ErosionAllM.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-ErosionAllM.DeltaAIC);
ErosionAllM=ErosionAllM(si,:);
ErosionAllM.Rank=si;
ErosionAllM.L=exp(0.5*ErosionAllM.DeltaAIC);
ErosionAllM.Weight=ErosionAllM.L./sum(ErosionAllM.L);

% MHI Net--------------------------------------------
y=(Data2.PercentNetChange(MHIa));

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

NetAllM=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s',p');
NetAllM.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-NetAllM.DeltaAIC);
NetAllM=NetAllM(si,:);
NetAllM.Rank=si;
NetAllM.L=exp(0.5*NetAllM.DeltaAIC);
NetAllM.Weight=NetAllM.L./sum(NetAllM.L);

%NWHI Accretion---------------------------------------
X=All(cumsum(index),:); 
X=X(NWHIa,:); %all environmental variables
y=log(Data2.AccRate(NWHIa));

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

AccretionAllN=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s',p');
AccretionAllN.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-AccretionAllN.DeltaAIC);
AccretionAllN=AccretionAllN(si,:);
AccretionAllN.Rank=si;
AccretionAllN.L=exp(0.5*AccretionAllN.DeltaAIC);
AccretionAllN.Weight=AccretionAllN.L./sum(AccretionAllN.L);

%NWHI EROSION-------------------------------------------------------
y=log(Data2.ErRate(NWHIa));

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

ErosionAllN=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s',p');
ErosionAllN.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-ErosionAllN.DeltaAIC);
ErosionAllN=ErosionAllN(si,:);
ErosionAllN.Rank=si;
ErosionAllN.L=exp(0.5*ErosionAllN.DeltaAIC);
ErosionAllN.Weight=ErosionAllN.L./sum(ErosionAllN.L);

%NWHI Net--------------------------------
y=Data2.PercentNetChange(NWHIa);

for i=1:width(X)
    mdl=fitlm(table2array(X(:,i)),y);
    loglik(i)=mdl.LogLikelihood;
    AIC(i)=mdl.ModelCriterion.AIC;
    R2(i)=mdl.Rsquared.Ordinary;
    Model{i}=names{i};
    s(i)=mdl.Coefficients.Estimate(2);
    p(i) = mdl.coefTest;
end

NetAllN=table(Model',loglik',AIC',(min(AIC)-AIC)',R2',s', p');
NetAllN.Properties.VariableNames={'Model','LogLik','AIC','DeltaAIC','R2','slope','p'};
[d,si]=sort(-NetAllN.DeltaAIC);
NetAllN=NetAllN(si,:);
NetAllN.Rank=si;
NetAllN.L=exp(0.5*NetAllN.DeltaAIC);
NetAllN.Weight=NetAllN.L./sum(NetAllN.L);

%% Bar plots of the Model Selection

%set the colors for physics, biology, and chemistry
cols=[repmat([1 0.4 0.6],[7,1]);repmat([0 1 1],[7,1]);repmat([1 1 0],[9,1])];
j=[3,2,1]; %needed to get the bar colors right 

figure;
ax1=subplot(3,3,1);
%HA accretion
xn=AccretionAll;
xr=xn.Weight;
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

for i=1:length(xn.slope) %show whether the slope was + or -
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold


xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
title('Hawaiian Archipelago')
ylabel('Secondary Accretion')
xlim([0 1])

ax2=subplot(3,3,2);
%MHI accretion
xn=AccretionAllM;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
title('MHI')
ax2.XLim=[0 1];

ax3=subplot(3,3,3);
%NWHI accretion
xn=AccretionAllN;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
title('NWHI')
xlim([0 1])

ax4=subplot(3,3,4);
%HA Erosion
xn=ErosionAll;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
ylabel('Bioerosion')
xlim([0 1])
ax5=subplot(3,3,5);
%MHI Erosion
xn=ErosionAllM;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
xlim([0 1])

ax6=subplot(3,3,6);
%NWHI Erosion
xn=ErosionAllN;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
xlim([0 1])
%All net
ax7=subplot(3,3,7);
xn=NetAll;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
ylabel('Net Change')
xlabel('AIC_w')
xlim([0 1])
%MHI Net
ax8=subplot(3,3,8);
xn=NetAllM;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold
%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
xlabel('AIC_w')
xlim([0 1])
%NWHI Net
ax9=subplot(3,3,9);
xn=NetAllN;
xr=xn.Weight;

for i=1:length(xn.slope)
if(sign(xn.slope(i))==-1); 
    u(i,1)='-'; 
else u(i,1)='+';
end
end

R=cellstr(num2str(round(xn.R2,2)));
R=strcat('\bf ',R,' (',num2str(u),')'); %make bold

%sort the colors so that they match
[t, ia, ib] = intersect(xn.Model,names,'stable');
c=cols(ib,:);

xr=flipud(xr(1:3));
hold on
for i = 1:numel(xr)
    barh(i, xr(i), 'facecolor', c(j(i),:));
end
set(gca,'XTick',0:0.25:1,'xticklabel',{0,[],0.5,[],1},'YTick', 1:3, 'YTickLabel', flipud(xn.Model(1:3)))
text(repmat(0.3,[3 1]),1:3, (flipud(R(1:3))))
xlabel('AIC_w')
xlim([0 1])




