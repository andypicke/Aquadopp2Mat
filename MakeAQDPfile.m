function aqdp=MakeAQDPfile(datadir,fnamebase,whbins,whsamp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function aqdp=MakeAQDPfile(datadir,fnamebase)
%
% Read Aquadopp data files and put all
% data into a matlab structure 'aqdp'
%
% AQDP files are in folder 'datadir' and are named
% 'fnamebase.hdr','fnamebase.v1' etc.
%
%
% Original AP 2 Aug 2012
% 16 Sept - Adding option to get only some bins (to reduce file size for
% large deployments or if only interested in 1 beam).
%
% 4 Oct 2012 AP - Also get amplitude & correlation data
% 4 Oct AP - add option to only get specified time range (whsamp)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

% check if all the files we need are there

if ~exist(fullfile(datadir,[fnamebase '.hdr']))
    error('header file doesnt exist')
elseif ~exist(fullfile(datadir,[fnamebase '.sen']))
    error('.sen file doesnt exist')
elseif ~exist(fullfile(datadir,[fnamebase '.v1']))
    error('.v1 file doesnt exist')
elseif ~exist(fullfile(datadir,[fnamebase '.v2']))
    error('.v2 file doesnt exist')
elseif ~exist(fullfile(datadir,[fnamebase '.v3']))
    error('.v3 file doesnt exist')
end

fid=fopen(fullfile(datadir,[fnamebase '.hdr']));

if fid==-1
    error('Cannot open header file')
end



%~~~~
%Find the line in hdr file for Coordinate system
disp('Reading header file for coordinate system ')
foundit=0;
whline=1;
%
while foundit==0
    %
    tline=fgetl(fid);
    %
    %
    if length(tline>0) % some lines are empty...
        %
        if strcmp(tline(1:10),'Coordinate') % find line that starts with 'Coordinate'
            foundit=1;
        else
        end
        %
    end
    %
    whline=whline+1;
    %
end

% Options will be 'BEAM','ENU',or 'XYZ'
if strcmp(tline(end-3:end),'BEAM')
    aqdp.CoordSys_orig=tline(end-3:end);
else
    aqdp.CoordSys_orig=tline(end-2:end);
end

disp(['Coordinate system is ' aqdp.CoordSys_orig ])

%%
%~~~~
% Get the serial number

foundit=0;
whline=1;
%
while foundit==0
    %
    tline=fgetl(fid);
    %
    %
    if length(tline>0) % some lines are empty...
        %
        if strcmp(tline(1:10),'Serial num') % find line that starts with
            foundit=1;
        else
        end
        %
    end
    %
    whline=whline+1;
    %
end

aqdp.SerialNum=str2num(tline(end-13:end));
%%
%~~~

%
%~~~~
% Get the instrument transform matrix (needed for coordinate
% system transform
disp('getting transform matrix')
foundit=0;
whline=1;
%
while foundit==0
    %
    tline=fgetl(fid);
    %
    %
    if length(tline>0) % some lines are empty...
        %
        if strcmp(tline(1:10),'Transforma') % find line that starts with
            foundit=1;
            disp('Found transform matrix')
        else
        end
        %
    end
    %
    whline=whline+1;
    %
end

aqdp.T=nan*ones(3,3);
aqdp.T(1,1:3)=str2num(tline(end-21:end));
tline=fgetl(fid);
aqdp.T(2,1:3)=str2num(tline(end-21:end));
tline=fgetl(fid);
aqdp.T(3,1:3)=str2num(tline(end-21:end));

%~~~
%


%
% Now load sensor data (time, heading, pitch, roll etc
%out=LoadAQDPdat(fullfile(datadir,[cruise '.sen']));
disp('loading sensor data (time, pitch, roll etc')
out=importdata(fullfile(datadir,[fnamebase '.sen']));

%
day=out(whsamp,2);
month=out(whsamp,1);
year=out(whsamp,3);
hour=out(whsamp,4);
min=out(whsamp,5);
sec=out(whsamp,6);

%
dnum=nan*ones(1,length(day));
yday=dnum;

% convert to datenum and yday
for jj=1:length(day)
    %    DisplayProgress(jj,1000)
    dnum(jj)=datenum(year(jj),month(jj),day(jj),hour(jj),min(jj),sec(jj));
end

aqdp.dtnum=dnum;
aqdp.yday=datenum2yday(dnum);

clear month day year hour min sec


% get heading,pitch,roll
disp('getting pitch and roll data')
aqdp.hdg=out(whsamp,13);
aqdp.pitch=out(whsamp,14);
aqdp.roll=out(whsamp,15);

% get pressure and temp.
disp('getting pressure and temp data')
aqdp.p=out(whsamp,16);
%dtsec=nanmean(diffs(aqdp.yday(:)))*86400;
%aqdp.dpdt=diffs(aqdp.p(:))./dtsec;
aqdp.t=out(whsamp,17);
clear out
%
% Now load velocity data
%
% Velocity columns are: burst / ensemble / bin1 / bin2 ....Nbins
%
clear out
disp('loading v1 ')
out=importdata(fullfile(datadir,[fnamebase '.v1']));
disp('v1 loaded; saving data ')

[mm,nn]=size(out);
mm=length(whsamp);

if ~exist('whbins','var')
    whbins=3:nn;
else
    whbins=whbins+2; %
end

%whbins
%pause

aqdp.burst=nan*ones(mm,1);
aqdp.ensemble=nan*ones(mm,1);

aqdp.burst=out(whsamp,1);
aqdp.ensemble=out(whsamp,2);


if strcmp(aqdp.CoordSys_orig,'BEAM')
    aqdp.v1=nan*ones(mm,nn);
    aqdp.v1=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'ENU')
    aqdp.u=nan*ones(mm,nn);
    aqdp.u=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'XYZ')
    aqdp.x=nan*ones(mm,nn);
    aqdp.x=out(whsamp,whbins);
end

clear out

disp('loading v2 ')
out=importdata(fullfile(datadir,[fnamebase '.v2']));
disp('v2 loaded; saving data ')
if strcmp(aqdp.CoordSys_orig,'BEAM')
    aqdp.v2=nan*ones(mm,nn);
    aqdp.v2=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'ENU')
    aqdp.v=nan*ones(mm,nn);
    aqdp.v=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'XYZ')
    aqdp.y=nan*ones(mm,nn);
    aqdp.y=out(whsamp,whbins);
end
%aqdp.v2=out(:,3:end);

clear out

disp('loading v3 ')
out=importdata(fullfile(datadir,[fnamebase '.v3']));
disp('v3 loaded; saving data ')
if strcmp(aqdp.CoordSys_orig,'BEAM')
    aqdp.v3=nan*ones(mm,nn);
    aqdp.v3=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'ENU')
    aqdp.w=nan*ones(mm,nn);
    aqdp.w=out(whsamp,whbins);
elseif strcmp(aqdp.CoordSys_orig,'XYZ')
    aqdp.z=nan*ones(mm,nn);
    aqdp.z=out(whsamp,whbins);
end
clear out

% amplitudes for beam 1
disp('loading a1 ')
out=importdata(fullfile(datadir,[fnamebase '.a1']));
aqdp.a1=nan*ones(mm,nn);
aqdp.a1=out(whsamp,whbins);
clear out

% amplitudes for beam 2
disp('loading a2 ')
out=importdata(fullfile(datadir,[fnamebase '.a2']));
aqdp.a2=nan*ones(mm,nn);
aqdp.a2=out(whsamp,whbins);
clear out

% amplitudes for beam 3
disp('loading a3 ')
out=importdata(fullfile(datadir,[fnamebase '.a3']));
aqdp.a3=nan*ones(mm,nn);
aqdp.a3=out(whsamp,whbins);
clear out
    
% correlations for beam 1
disp('loading c1 ')
out=importdata(fullfile(datadir,[fnamebase '.c1']));
aqdp.c1=nan*ones(mm,nn);
aqdp.c1=out(whsamp,whbins);
clear out

% correlations for beam 1
disp('loading c2 ')
out=importdata(fullfile(datadir,[fnamebase '.c2']));
aqdp.c2=nan*ones(mm,nn);
aqdp.c2=out(whsamp,whbins);
clear out

% correlations for beam 1
disp('loading c3 ')
out=importdata(fullfile(datadir,[fnamebase '.c3']));
aqdp.c3=nan*ones(mm,nn);
aqdp.c3=out(whsamp,whbins);
clear out

aqdp.processed=[date ' with ' mfilename];

disp('Done!')

return
%%
