function [geostrophic_velocity, mid_lat, mid_long] = gsw_geostrophic_velocity(geo_str,long,lat,p)

% USAGE:  
%  [geostrophic_velocity, mid_lat, mid_long] = gsw_geostrophic_velocity(geo_strf,long,lat,p)
%
% DESCRIPTION:
%  Calculates geostrophic velocity relative to the sea surface, given a
%  geostrophic streamfunction and the position (longitude, latitude and
%  pressure (long, lat & p)) of each station in sequence along an ocean
%  section.  The data can be from a single isobaric or "density" surface,
%  or from a series of such surfaces. 
%
% INPUT:
%  geo_strf = geostrophic streamfunction. This geostrophic streamfunction
%             can be any of, for example,
%             (1) geo_strf_dyn_height (in an isobaric surface)
%             (2) geo_strf_Montgomery (in a specific volume anomaly surface)
%             (3) geo_strf_Cunninhgam (in an approximately neutral surface
%                 e.g. a potential denisty surface).
%             (4) geo_strf_McD_Klocker (in an approximately neutral surface
%                 e.g. a potential denisty surface, a Neutral Density
%                 surface or an omega surface (Klocker et al., 2009)).
%
%  long    =  Longitude in decimal degrees                   [ 0 ... +360 ]
%                                                      or [ -180 ... +180 ]
%  lat     =  Latitude in decimal degrees north             [ -90 ... +90 ]
%  
%  There needs to be more than one station.  The input geo_strf has
%  dimensions (M("bottles") x N(stations)); that is, geo_strf has 
%  dimensions (M(surfaces)  x N(stations)). 
%
%  Note. The ith "bottle" of each station (ie. the ith row of geo_strf) 
%    must be on the same ith surface, whether that surface be, 
%     (1) an isobaric surface, 
%     (2) a specific volume anomaly surface, 
%         or some type of approximately neutral surface (cases (3) & (4)). 
%
%  lat & long need to have dimensions 1xN or MxN, where geo_strf is MxN.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where geo_strf is MxN.
%
% OPTIONAL INPUT:
%  p     =  sea pressure ( default is 0 )                          [ dbar ]
%           (i.e. absolute pressure - 10.1325 dbar)
%  Note. This optional input is used to obtain an accurate distance,
%    "dist", taking into account that the radius from the centre of the 
%    Earth depends on the depth below the sea surface.               
%
%
% OUTPUT:
%  geostrophic_velocity = geostrophic velocity RELATIVE to the sea surface.
%                         It has dimensions (Mx(N-1)) and the relative 
%                         geostrophic velocity is located at mid_long, 
%                         mid_lat (and at the mid-point pressure).  
%
%  mid_long             = mid point longitude
%                         (the range corresponds to that entered)
%
%  mid_lat              = mid point latitude,               [ -90 ... +90 ]
%                         (in decimal degrees North)
%
%
% AUTHOR: Paul Barker, Trevor McDougall and Phil Morgan 
%                                                     [ help_gsw@csiro.au ]
% VERSION NUMBER: 2.0 (14th September, 2010)
%
% REFERENCES:
%  Cunningham, S. A., 2000: Circulation and volume flux of the North 
%   Atlantic using syoptic hydrographic data in a Bernoulli inverse.
%   J. Marine Res., 58, 1-35. 
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See sections 3.27 - 3.3.30 of this TEOS-10 Manual. 
%
%  Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable 
%   for the world’s oceans. Journal of Physical Oceanography, 27, 237-263.
%
%  Klocker, A., T. J. McDougall and D. R. Jackett, 2009: A new method for
%   forming approximately neutral surfaces.  Ocean Sci., 5, 155-172. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  McDougall, T. J. and A. Klocker, 2010: An approximate geostrophic 
%   streamfunction for use in density surfaces.  Ocean Modelling, 32, 
%   105-117. 
%    See Eqn. (62), of this paper, for definition of the McDougall-Klocker
%    geostrophic streamfunction. 
%
%  Montgomery, R. B., 1937: A suggested method for representing gradient 
%   flow in isentropic surfaces.  Bull. Amer. Meteor. Soc. 18, 210-212.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 | nargin == 4)
   error('gsw_geostrophic_velocity:  Requires three or four inputs')
end %if

[mg,ng] = size(geo_str);
[mlo,nlo] = size(long);
[mla,nla] = size(lat);

if (nlo == 1) & (mlo == 1)             
    error('gsw_geostrophic_velocity: need more than one station')
elseif (ng == nlo) & (mlo == 1)                       % long is row vector,
    long = long(ones(1,mg), :);                    % copy down each column.
elseif (ng == mlo) & (nlo == 1)          % long is a transposed row vector,
    long = long';                                         % transposed then
    long = long(ones(1,mg), :);                    % copy down each column.
elseif (mg == mlo) & (ng == nlo)
    % ok
else
    error('gsw_geostrophic_velocity: Inputs array dimensions arguments do not agree')
end %if

[Iwest] = find(long < 0);
if ~isempty(Iwest)
    long(Iwest) = long(Iwest) + 360;
end

if (nla == 1) & (mla == 1)            
    error('gsw_geostrophic_velocity: need more than one station')
elseif (ng == nla) & (mla == 1)                        % lat is row vector,
    lat = lat(ones(1,mg), :);                      % copy down each column.
elseif (ng == mla) & (nla == 1)          % long is a transposed row vector,
    lat = lat';                                           % transposed then
    lat = lat(ones(1,mg), :);                      % copy down each column.
elseif (mg == mla) & (ng == nla)
    % ok
else
    error('gsw_geostrophic_velocity: Inputs array dimensions arguments do not agree')
end %if

if exist('p','var')
    [mp,np] = size(p);
    if (mp == 1) & (np == 1)           % p scalar - fill to size of geo_str
        p = p*ones(size(geo_str));
    elseif (ng == np) & (mp == 1)                        % p is row vector,
        p = p(ones(1,mg), :);                      % copy down each column.
    elseif (mg == mp) & (np == 1)                     % p is column vector,
        p = p(:,ones(1,ng));                        % copy across each row.
    elseif (ng == mp) & (np == 1)           % p is a transposed row vector,
        p = p';                                           % transposed then
        p = p(ones(1,mg), :);                      % copy down each column.
    elseif (mg == mp) & (ng == np)
        % ok
    else
        error('gsw_geostrophic_velocity: Inputs array dimensions arguments do not agree')
    end %if
else
    p = zeros(size(geo_str));
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

dist = gsw_distance(long,lat,p);          % Note that dist is in m (not km)

mid_lat = 0.5*(lat(:,1:ng-1) + lat(:,2:ng));

f = gsw_f(mid_lat); % This gets the Coriolis parameter

geostrophic_velocity = (geo_str(:,2:ng) - geo_str(:,1:ng-1))./(dist.*f);

mid_long = 0.5*(long(:,1:ng-1) + long(:,2:ng));
diff_long = (long(1,1:ng-1) - long(1,2:ng));

long2 = long;
[Iwest2] = find(long2 > 180);
long2(Iwest2) = long2(Iwest2) - 360;
mid_long2 = 0.5*(long2(:,1:ng-1) + long2(:,2:ng));
diff_long2 = (long2(1,1:ng-1) - long2(1,2:ng));

[dummy,gmc] = min([abs(diff_long); abs(diff_long2)]);

[Igmc] = find(gmc == 2);
if ~isempty(Igmc)
    mid_long(Igmc) = mid_long2(Igmc);
end

[Iwest2] = find(mid_long < 0);
if ~isempty(Iwest2)
    mid_long(Iwest2) = mid_long(Iwest2) + 360;
end

if ~isempty(Iwest)
    [Iwestrestore] = find(mid_long > 180);
    mid_long(Iwestrestore) = mid_long(Iwestrestore) - 360;
end

% Note. This geostrophic velocity difference, v - v0, when positive, is 
%   directed to the left of the horizontal vector which points from one
%   cast to the next.  

end

