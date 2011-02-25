function [geo_strf_McD_Klocker, in_funnel] = gsw_geo_strf_McD_Klocker(SA,CT,p,Neutral_Density,p_Neutral_Density,A)

% gsw_geo_strf_McD_Klocker      McDougall-Klocker geostrophic streamfunction
%==========================================================================
%
% USAGE:
%  [geo_strf_McD_Klocker, in_funnel] = gsw_geo_strf_McD_Klocker(SA,CT,p,Neutral_Density,p_Neutral_Density,A)
%
% DESCRIPTION:
%  Calculates the McDougall-Klocker geostrophic streamfunction (see Eqn. 
%  (3.30.1) of IOC et al. (2010)).  This is the geostrophic streamfunction 
%  for the difference between the horizontal velocity at the pressure 
%  concerned, p, and the horizontal velocity at the sea surface.  It is 
%  designed to be used as the geostrophic streamfunction in an 
%  approximately neutral surface (such as a Neutral Density surface, a 
%  potential density surface or an omega surface (Klocker et al. (2009)).  
%  Reference values of Absolute Salinity, Conservative Temperature and 
%  pressure are found by interpolation of a one-dimensional look-up table, 
%  with the interpolating variable being Neutral Density (gamma_n).  This
%  function calculates specific volume anomaly using the computationally 
%  efficient 25-term expression for specific volume in terms of SA, CT and
%  p (McDougall et al., 2010). 
%
%  The first three input arguments are a series of vertical profiles, while
%  the last two argumnents pertain to the (usually relatively few) surfaces
%  on which the McDougall-Klocker geostrophic streamfunction is to be 
%  calculated.  These last two input arguments, Neutral_Density and 
%  p_Neutral_Density, are the Neutral density label and the pressure of 
%  each of the (usually relately few) surfaces.  p_Neutral_Density is the 
%  series of pressures where the surfaces intersect the vertical profiles. 
%  These surfaces do not have to be the very best approximately neutral 
%  surfaces; rather the onus is on the user to use a surface that is 
%  sufficiently neutral for his/her purpose.   The input variable 
%  "Neutral_Density" is used to find reference values of SA, CT and p 
%  by vertcal interpolation down a single reference cast.  As an 
%  alternative to the user supplying Neutral Density for this purpose, 
%  the code allows for sigma_2 to be used as the vertical interpolating
%  variable instead of Neutral Density.     
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  CT       =  Conservative Temperature                           [ deg C ]
%  p        =  sea pressure                                        [ dbar ]
%              (ie. absolute pressure - 10.1325 dbar)
%
%  Neutral_Density   =  Neutral Density anomaly                  [ kg/m^3 ]
%              (ie. Neutral Density minus 1000 kg/m^3)
%  p_Neutral_Density =  pressure of the Neutral_Density surface.
%
%  A        =  if nothing is entered the programme defaults to "Neutral
%              Density" as the vertical interpolating variable. 
%           = 's2' or 'sigma2', for sigma_2 as the vertical interpolating
%              variable. 
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  Neutral_Density & p_Neutral_Density need to have the same dimensions,
%  and they need to have dimensions BxN, where B is the number of surfaces.
%
% OUTPUT:
%  geo_strf_McD_Klocker = McDougall & Klocker (2010)             [ m^2/s^2 ]
%                        geostrophic streamfunction
%
%  in_funnel  =  0, if SA, CT and p are outside the "funnel" 
%             =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR:  Trevor McDougall and Paul Barker
%
% VERSION NUMBER: 2.0 (1st September, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.30 of this TEOS-10 Manual.
%
%  Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable
%   for the world’s oceans. Journal of Physical Oceanography, 27, 237-263.
%
%  Klocker, A., T. J. McDougall and D. R. Jackett, 2009: A new method 
%   for forming approximately neutral surfaces.  Ocean Sci., 5, 155-172. 
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
%    The McDougall-Klocker geostrophic streamfunction is defined in
%   Eqn. (62) of this paper.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5 | nargin == 6 )
    error('gsw_geo_strf_McD_Klocker:  Requires five or six inputs')
end %if

if ~exist('A','var')
    A = 'gn';
elseif ~ischar(A)
    A = 'gn';
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_McD_Klocker: SA & CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    error('gsw_geo_strf_McD_Klocker: needs more than one pressure');
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_geo_strf_McD_Klocker: Inputs array dimensions arguments do not agree')
end %if

[mgn,ngn] = size(Neutral_Density);
[mpgn,npgn] = size(p_Neutral_Density);

if (mgn~=mpgn) | (ngn~=npgn)
    error('gsw_geo_strf_McD_Klocker: Neutral_Density & p_Neutral_Density need to have the same dimensions')
end

if mgn == 1 & ms == 1
    Neutral_Density = Neutral_Density(:);
    p_Neutral_Density = p_Neutral_Density(:);
end

if ms == 1  % row vector
    p  =  p(:);
    CT  =  CT(:);
    SA  =  SA(:);
    Neutral_Density = Neutral_Density(:);
    p_Neutral_Density = p_Neutral_Density(:);
    transposed = 1;
else
    transposed = 0;
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mgn,ngn] = size(Neutral_Density);
[mpgn,npgn] = size(p_Neutral_Density);

if ngn ~= ns
    error('gsw_geo_strf_McD_Klocker: SA & Neutral_Density need to have the same number of profiles')
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;
cp0 = 3991.86795711963;           % from Eqn. (3.3.3) of IOC et al. (2010).

SA_iref_cast = nan(size(Neutral_Density));
CT_iref_cast = nan(size(Neutral_Density));
p_iref_cast = nan(size(Neutral_Density));

[Inn] = find(~isnan(Neutral_Density));
[SA_iref_cast(Inn),CT_iref_cast(Inn),p_iref_cast(Inn)] = gsw_interp_McD_Klocker(Neutral_Density(Inn),A);

dummy = cat(1,p,p_Neutral_Density); % this combines the profile pressures with the neutral density pressures
p_plusnd = sort(dummy);             % this sorts the pressures down the profile into decending order

dyn_height_nd = nan(mpgn,ns);
SA_nd = nan(mpgn,ns);
CT_nd = nan(mpgn,ns);
in_funnel= nan(mpgn,ns);

for Iprofile = 1:ns
    [Inn] = ~isnan(p_plusnd(:,Iprofile));
    p_plus = unique(p_plusnd(Inn,Iprofile));
    [InnSACT] = ~isnan(SA(:,Iprofile)) & ~isnan(CT(:,Iprofile));
    [SA_plus, CT_plus] = gsw_interp_SA_CT(SA(InnSACT,Iprofile),CT(InnSACT,Iprofile),p(InnSACT,Iprofile),p_plus);
    [dyn_height,in_funnel_plus] = gsw_geo_strf_dyn_height(SA_plus(:),CT_plus(:),p_plus);
    [dummy, Idata] = intersect(p_plus,p_Neutral_Density(:,Iprofile));
    dyn_height_nd(1:mpgn,Iprofile) = dyn_height(Idata);
    SA_nd(1:mpgn,Iprofile) = SA_plus(Idata);
    CT_nd(1:mpgn,Iprofile) = CT_plus(Idata);
    in_funnel(1:mpgn,Iprofile) = in_funnel_plus(Idata);
end

[Isurface] = find(p_Neutral_Density == 0);
p_Neutral_Density(Isurface) = NaN;

part1 = 0.5*db2Pa*(p_Neutral_Density -p_iref_cast).*(gsw_specvol_CT25(SA_nd,CT_nd,p_Neutral_Density) - ...
                                   gsw_specvol_CT25(SA_iref_cast,CT_iref_cast,p_Neutral_Density));

part2 = -0.225e-15*db2Pa*db2Pa*(CT_nd-CT_iref_cast).*(p_Neutral_Density-p_iref_cast).*(p_Neutral_Density-p_iref_cast);

part3 = dyn_height_nd - gsw_enthalpy_SSO_0_CT25(p_Neutral_Density) + ...
        gsw_enthalpy_CT25(SA_iref_cast,CT_iref_cast,p_Neutral_Density) - cp0*CT_iref_cast;

%--------------------------------------------------------------------------
% This function calculates the McDougall-Klocker streamfunction using the 
% computationally efficient 25-term expression for density in terms of SA, 
% CT and p.  If one wanted to compute this with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable 
% this.  Note that dynamic height will also need to be evaluated using the
% full Gibbs function.
% 
% part1 = 0.5*db2Pa*(p_Neutral_Density -p_iref_cast).*(gsw_specvol_CT(SA_nd,CT_nd,p_Neutral_Density) - ...
%                                    gsw_specvol_CT(SA_iref_cast,CT_iref_cast,p_Neutral_Density));
% part2 = -0.225e-15*db2Pa*db2Pa*(CT_nd-CT_iref_cast).*(p_Neutral_Density-p_iref_cast).*(p_Neutral_Density-p_iref_cast);
% SA_SO = 35.16504*ones(size(SA));
% CT_0 = zeros(size(CT));
% part3 = dyn_height_nd - gsw_enthalpy_CT(SA_SO,CT_0,p_Neutral_Density) + ...
%         gsw_enthalpy_CT(SA_iref_cast,CT_iref_cast,p_Neutral_Density) - cp0*CT_iref_cast;
% 
%---------------This is the end of the alternative code--------------------

geo_strf_McD_Klocker = part1 + part2 + part3;

if transposed
    geo_strf_McD_Klocker= geo_strf_McD_Klocker';
    in_funnel = in_funnel';
end %if

end


