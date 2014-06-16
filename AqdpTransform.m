function [aqdp] = AqdpTransform(aqdp,transform)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Do coordinate transformation for Aquadopp HR velocities.
% Assuming aqdp structure made from MakeAQDPfile.m
%
% Options for transform:
% 'be' beam to earth (ENU)
% 'bx' beam to XYZ (instrument frame)
% 'xe' XYZ to earth
%  etc...
%
% Naming convention:
% Beam velocity: v1,v2,v3
% XYZ velocity :  vx,vy,vz
% ENU velocity :  u,v,w
%
%
% Transformed velocities are added to the aqdp structure.
%
% The Xform is done in order: Roll,Pitch,heading
%
%
% ***
% Note velocities may not be accurate for pitch or roll >10 deg
% - At higher angles, the order of rotations and angle conventions
% (ie Euler angles) become important ...
% ***
%
% Modifed from ParadoppXform.m
% AP 3 Aug 2012
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Old documentation:
% ParadoppXform.m is a Matlab function adpated from Transform.m (available
% on the NortekAS forums) that transforms velocity data between
% beam, XYZ, and ENU coordinates.
%
% transform specifies the desired transformation and should be
% 'bx' (beam to XYZ), 'xb' (XYZ to beam), 'xe' (XYZ to ENU). Case doesn't
% matter
%
% beam coordinates are defined as the velocity measured along the three
% beams of the instrument.
%
% XYZ coordinates are defined relative to the instrument, consult the
% specific manual for this.
%
%
%
% ENU coordinates are defined in an earth coordinate system, where
% E represents the East-West component, N represents the North-South
% component and U represents the Up-Down component.
%
% Note that the transformation matrix must be recalculated every time
% the orientation, heading, pitch or roll changes.

% Transformation matrix for beam to xyz coordinates,
% the values can be found from the header file that is generated in the
% conversion of data to ASCII format


%T = PD.system.transform.T;   % Scale the transformation matrix correctly to floating point numbers (this is done already in some files...)

T=aqdp.T;

T_orig = T;

% If instrument is pointing down (bit 0 in status equal to 1)
% rows 2 and 3 must change sign
% NOTE: For the Vector the instrument is defined to be in
%       the UP orientation when the communication cable
%       end of the canister is on the top of the canister, ie when
%       the probe is pointing down.
%       For the other instruments that are used vertically, the UP
%       orientation is defined to be when the head is on the upper
%       side of the canister.
% we'll assume the system




%** don't need this?
% if isfield(PD,'setup.orientation')
%     orientation = PD.setup.orientation;
% else
%     orientation = 'd';
% end
% 
% if strcmpi(orientation,'d')
%     T([2 3],:) = -1 * T([2 3],:);
% end


Nsamples=length(aqdp.dtnum);
Nbins=size(aqdp.v1,2);

% First do beam<-->xyz transformations. For these, the rotation matrix
% is just T and is the same for all times.  S


% if ndims(Uin) == 2 % this is Vector data, much easier
%     if size(Uin,1) ~= 3
%         Uin = Uin';
%     end

switch lower(transform)
    case 'bx'
        % beam to XYZ coordinates
        %size(Uin2)
        %size(inv(T_orig))
        
        %
        disp('Doing BEAM --> XYZ transformation')
        
        % sample # X bin# 
        U1out=nan*ones(Nsamples,Nbins);
        U2out=nan*ones(Nsamples,Nbins);
        U3out=nan*ones(Nsamples,Nbins);

        for whs=1:Nsamples
            clear Uin Uout2
            Uin=[aqdp.v1(whs,:) ; aqdp.v2(whs,:) ; aqdp.v3(whs,:)];
            Uout2 = dothetransform(Uin,T_orig);
            U1out(whs,:)=Uout2(1,:);
            U2out(whs,:)=Uout2(2,:);
            U3out(whs,:)=Uout2(3,:);
        end
        
        aqdp.vx=U1out;
        aqdp.vy=U2out;
        aqdp.vz=U3out;
        
        clear U1out U2out U3out
        
        
    case 'xb'
        % XYZ to beam coordinates
                %
        disp('Doing XYZ --> BEAM transformation')
        
        if isfield(aqdp,'v1')
            disp('Warning: beam velocity fields already exist?')
        end
        
        % sample # X bin# 
        U1out=nan*ones(Nsamples,Nbins);
        U2out=nan*ones(Nsamples,Nbins);
        U3out=nan*ones(Nsamples,Nbins);
         
        for whs=1:Nsamples
            clear Uin Uout2
            Uin=[aqdp.v1(whs,:) ; aqdp.v2(whs,:) ; aqdp.v3(whs,:)];
            Uout2 = dothetransform(Uin,inv(T_orig));
            U1out(whs,:)=Uout2(1,:);
            U2out(whs,:)=Uout2(2,:);
            U3out(whs,:)=Uout2(3,:);
        end
        
        aqdp.v1=U1out;
        aqdp.v2=U2out;
        aqdp.v3=U3out;
        
        clear U1out U2out U3out
        
        
    case 'xe'
        
        disp('Doing XYZ --> ENU transformation')
        % XYZ to ENU coordinates - need to make new rotation matrix
        % for each time/sample.

        % sample # X bin# 
        U1out=nan*ones(Nsamples,Nbins);
        U2out=nan*ones(Nsamples,Nbins);
        U3out=nan*ones(Nsamples,Nbins);
         
        if ~isfield(aqdp,'vx')
            error('No XYZ fields. Need to convert beam to XYZ first (bx) ')
        end
        
        for whs=1:Nsamples
            clear Uin Uout2 hh pp rr H P R
            Uin=[aqdp.vx(whs,:) ; aqdp.vy(whs,:) ; aqdp.vz(whs,:)];

            % heading, pitch, and roll for this sample
                hh = pi*(aqdp.hdg(whs)-90)/180;
                pp = pi*aqdp.pitch(whs)/180;
                rr = pi*aqdp.roll(whs)/180;
                
                % Make heading matrix
                H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                
                % Make tilt matrix
                P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                    0             cos(rr)          -sin(rr);  ...
                    sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                
                % Make resulting transformation matrix
                R = H*P;
        
                Uout2= dothetransform(Uin,R);
            
            U1out(whs,:)=Uout2(1,:);
            U2out(whs,:)=Uout2(2,:);
            U3out(whs,:)=Uout2(3,:);
        end
  
        aqdp.u=U1out;
        aqdp.v=U2out;
        aqdp.w=U3out;
        
        clear U1out U2out U3out
        
        

%     case 'ex' %** Should just be inverse of xe transformation?
%         disp('Doing ENU-->XYZ transformation')
%         % heading, pitch and roll are the angles output in the data in
%         % degrees, recorded at 1 Hz, so for each N = PD.fs samples
%         % compute a new matrix
%         for i = 1:floor(length(Uin)/PD.fs)
%             if isnan(PD.ens(i*PD.fs)) ~= 1
%                 
%                 hh = pi*(PD.attitude.heading(i)-90)/180;
%                 pp = pi*PD.attitude.pitch(i)/180;
%                 rr = pi*PD.attitude.roll(i)/180;
%                 % Make heading matrix
%                 H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
%                 
%                 % Make tilt matrix
%                 P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
%                     0             cos(rr)          -sin(rr);  ...
%                     sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
%                 
%                 % Make resulting transformation matrix
%                 R = H*P*T_orig;
%                 
%                 % Given ENU velocities, xyz coordinates are calculated as
%                 Uout2(:,(i-1)*PD.fs+1:i*PD.fs) = dothetransform(Uin(:,(i-1)*PD.fs+1:i*PD.fs),T_orig*inv(R));
%                 %xyz = T_orig*inv(R)*enu;
%             else
%                 Uout2(:,(i-1)*PD.fs+1:i*PD.fs) = NaN(size(Uin,1),PD.fs);
%             end
%             
%         end
        
        
        %~~~~~~~~~~~~~~~~~
        
    case 'be'
        % beam to ENU coordinates
        
          disp('Doing BEAM-->ENU transformation')
        % heading, pitch and roll are the angles output in the data in
        % degrees, recorded at 1 Hz, so for each N = PD.fs samples
        % compute a new matrix
        
                % sample # X bin# 
        U1out=nan*ones(Nsamples,Nbins);
        U2out=nan*ones(Nsamples,Nbins);
        U3out=nan*ones(Nsamples,Nbins);
         
        if isfield(aqdp,'u')
            warning('ENU fields already exist in aqdp?')
        end
        
        if ~isfield(aqdp,'v1')
            error('No BEAM fields. ')
        end
        
        for whs=1:Nsamples
            DisplayProgress(whs,10000)
            clear Uin Uout2 hh pp rr H P R
            Uin=[aqdp.v1(whs,:) ; aqdp.v2(whs,:) ; aqdp.v3(whs,:)];

            % heading, pitch, and roll for this sample
                hh = pi*(aqdp.hdg(whs)-90)/180;
                pp = pi*aqdp.pitch(whs)/180;
                rr = pi*aqdp.roll(whs)/180;
                
                % Make heading matrix
                H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                
                % Make tilt matrix
                P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                    0             cos(rr)          -sin(rr);  ...
                    sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                
                % Make resulting transformation matrix
                R = H*P*T_orig;
%                Uout2(:,(i-1)*PD.fs+1:i*PD.fs) = dothetransform(Uin(:,(i-1)*PD.fs+1:i*PD.fs),R);
                 Uout2= dothetransform(Uin,R);
                 
            U1out(whs,:)=Uout2(1,:);
            U2out(whs,:)=Uout2(2,:);
            U3out(whs,:)=Uout2(3,:);
        end
  
        aqdp.u=U1out;
        aqdp.v=U2out;
        aqdp.w=U3out;
        
        clear U1out U2out U3out        

            
      % Earth to beam. this is just inverse of beam to earth...
    case 'eb'
        % heading, pitch and roll are the angles output in the data in
        % degrees, recorded at 1 Hz, so for each N = PD.fs samples
        % compute a new matrix
                % Given ENU velocities, beam coordinates are calculated as

                                % sample # X bin# 
        U1out=nan*ones(Nsamples,Nbins);
        U2out=nan*ones(Nsamples,Nbins);
        U3out=nan*ones(Nsamples,Nbins);
         
        if isfield(aqdp,'v1')
            warning('BEAM fields already exist in aqdp?')
        end
        
        if ~isfield(aqdp,'u')
            error('No ENU fields to transform from. ')
        end
        
        for whs=1:Nsamples
            clear Uin Uout2 hh pp rr H P R
            Uin=[aqdp.u(whs,:) ; aqdp.v(whs,:) ; aqdp.w(whs,:)];

            % heading, pitch, and roll for this sample
                hh = pi*(aqdp.hdg(whs)-90)/180;
                pp = pi*aqdp.pitch(whs)/180;
                rr = pi*aqdp.roll(whs)/180;
                
                % Make heading matrix
                H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1];
                
                % Make tilt matrix
                P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
                    0             cos(rr)          -sin(rr);  ...
                    sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)];
                
                % Make resulting transformation matrix
                R = H*P*T_orig;
                Uout2= dothetransform( Uin,inv(R) );
                 
            U1out(whs,:)=Uout2(1,:);
            U2out(whs,:)=Uout2(2,:);
            U3out(whs,:)=Uout2(3,:);
        end
  
        aqdp.v1=U1out;
        aqdp.v2=U2out;
        aqdp.v3=U3out;
        
        clear U1out U2out U3out        
              
end


function [Vout] = dothetransform(Vin,T)
Vout = T*Vin;

% % heading, pitch and roll are the angles output in the data in degrees
% hh = pi*(PD.transform.heading-90)/180;
% pp = pi*PD.transform.pitch/180;
% rr = pi*PD.transform.roll/180;

% % Make heading matrix
% H = [cos(hh) sin(hh) 0; -sin(hh) cos(hh) 0; 0 0 1]

% % Make tilt matrix
% P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp);...
%       0             cos(rr)          -sin(rr);  ...
%       sin(pp) sin(rr)*cos(pp)  cos(pp)*cos(rr)]

% % Make resulting transformation matrix
% R = H*P*T;


% beam is beam coordinates, for example beam = [0.23 ; -0.52 ; 0.12]
% enu is ENU coordinates
%{
% Given beam velocities, ENU coordinates are calculated as
enu = R*beam

% Given ENU velocities, beam coordinates are calculated as
beam = inv(R)*enu


% Transformation between beam and xyz coordinates are done using
% the original T matrix
xyz = T_orig*beam
beam = inv(T_orig)*xyz

% Given ENU velocities, xyz coordinates are calculated as
xyz = T_orig*inv(R)*enu
%}