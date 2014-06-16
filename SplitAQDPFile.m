function [Nprofiles,outdir]=SplitAQDPFile(aqdp,outdir,fnamebase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function SplitAQDPFile(infile)
%
% Splits Aquadopp data file (made with MakeAQDPfile.m)
% into one strucutre/file for each up or down profile.
%
% *Originally written to integrate Aquadopp data with MP data.
% Now instead use SplitAqdpMP.m , which uses the MP profile
% times to ensure that MP and aquadopp profiles are matched
% up correctly.*
%
% Uses aquadopp dp/dt to determine start and end
% of up/down profiles.  Data for each profile is saved to
% a structure 'a'
%
% Saves individual files to folder 'outdir'
% Out file format is 'aqdp_cruise_01.mat' etc
%
% This might bog down if the MP has trouble crawling and 
% dpdt is funky...

%
% AP 2 Aug 2012
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%%

%fnamebase=[cruise '_' station '_' MPsn]

saveoutfiles=1

dpdt_thresh=0.1
dpdt=diffs(aqdp.p);
dpdt=smooth(dpdt,35);

figure(33);clf
plot(aqdp.yday,dpdt)

% use dpdt to find start of up and down profiles
iddn=find(dpdt>dpdt_thresh);
idup=find(dpdt<-dpdt_thresh);

bdn=FindContigSeq(iddn);
bup=FindContigSeq(idup);


whprf=1
ind=1
inddn=1
indup=1

pmax=nanmax([bdn.N bup.N]);

while whprf<pmax   
clear a ydayst_up ydayst_dn 

disp(['working on profile # ' num2str(whprf)])
    
ydayst_up=aqdp.yday(bup.first(indup));
ydayst_dn=aqdp.yday(bdn.first(inddn));

clear idprf

if ydayst_dn < ydayst_up
    idprf=bdn.first(inddn) : bdn.last(inddn);
    inddn=inddn+1;
else
    idprf=bup.first(indup) : bup.last(indup);
    indup=indup+1;
end

% exclude any very short profiles caused by spikes in dpdt
if idprf>5

a.prfnum=whprf;
a.ydaystart=aqdp.yday(idprf(1));
a.ydayend=aqdp.yday(idprf(end));
a.dtnum=aqdp.dtnum(idprf);
a.yday=aqdp.yday(idprf);
a.p=aqdp.p(idprf);

% figure(1);clf
% plot(aqdp.yday,aqdp.p,'k')
% hold on
% plot(a.yday,a.p,'ro-')
% title(['profile ' num2str(inddn)])
%pause


dt=nanmean(diffs(a.yday))*86400;
a.dpdt=aqdp.dpdt(idprf)./dt;

if nanmean(diffs(a.p))>0
    a.whdir='dn';
elseif nanmean(diffs(a.p))<0
    a.whdir='up';
end
    

if isfield(aqdp,'v1')
a.v1=aqdp.v1(idprf,:);
a.v2=aqdp.v2(idprf,:);
a.v3=aqdp.v3(idprf,:);
end

if isfield(aqdp,'u')
a.u=aqdp.u(idprf,:);
a.v=aqdp.v(idprf,:);
a.w=aqdp.w(idprf,:);
end


a.heading=aqdp.hdg(idprf);
a.pitch=aqdp.pitch(idprf);
a.roll=aqdp.roll(idprf);

a.CoordSys_orig=aqdp.CoordSys_orig;

a.Made=[date ' w/ ' mfilename ];
%a.cruise=cruise;
%a.station=station;
%a.MPsn=MPsn;

% save file

if saveoutfiles==1
%outdir=fullfile(basedir,'out');
%outdir='/Users/Andy/Cruises_Research/AquaDoppTests/outprofiles/'

if ~isdir(outdir)
disp('Out directory does not exist, creating now ')
mkdir(outdir)
else
end
    
fname=['aqdp_' fnamebase '_' sprintf('%03.0f',whprf) '.mat']
%ffull=fullfile(outdir,fname);
save(fullfile(outdir,fname),'a')
end
%pause

% increment counters
whprf=whprf+1;
%ind=ind+1;

a;

end

end

Nprofiles=whprf-1;

return

%%
