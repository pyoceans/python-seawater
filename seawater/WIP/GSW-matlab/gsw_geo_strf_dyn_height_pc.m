function [geo_strf_dyn_height_pc, p_mid, in_funnel] = gsw_geo_strf_dyn_height_pc(SA,CT,delta_p)

% gsw_geo_strf_dyn_height_pc                     dynamic height anomaly for
%                                               piecewise constant profiles
%==========================================================================
%
% USAGE:  
%  [geo_strf_dyn_height_pc, p_mid, in_funnel] = gsw_geo_strf_dyn_height_pc(SA,CT,delta_p)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume 
%  anomaly from the the sea surface pressure (0 Pa) to the pressure p.
%  This function, gsw_geo_strf_dyn_height_pc, is to used when the Absolute
%  Salinity and Conservative Temperature are piecewise constant in the
%  vertical over sucessive pressure intervals of delta_p (such as in
%  a forward "z-coordinate" ocean model).  "geo_strf_dyn_height_pc" is the 
%  dynamic height anomaly with respect to the sea surface.  That is, 
%  "geo_strf_dyn_height_pc" is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned, p,
%  and the horizontal velocity at the sea surface.  Dynamic height anomaly 
%  is the geostrophic streamfunction in an isobaric surface.  The reference
%  values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
%  and CT = 0 deg C.  The output values of geo_strf_dyn_height_pc are given 
%  at the mid-point pressures, p_mid, of each layer in which SA and CT are 
%  vertically piecewice constant(pc).  This function calculates specific 
%  volume anomaly using the computationally efficient 25-term expression 
%  for specific volume of McDougall et al. (2010). 
%
% INPUT:
%  SA       =  Absolute Salinity                                   [ g/kg ]
%  CT       =  Conservative Temperature                           [ deg C ]
%  delta_p  =  difference in sea pressure between the deep and     [ dbar ]
%              shallow extents of each layer in which SA and CT
%              are vertically constant delta_p must be positive.
%              
%  Note. Sea pressure is absolute pressure minus 10.1325 dbar.
%
%  SA & CT need to have the same dimensions.
%  delta_p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  geo_strf_dyn_height_pc  =  dynamic height anomaly            [ m^2/s^2 ]
%  p_mid                   =  mid-point pressure in each layer     [ dbar ]
%
%  in_funnel               =  0, if SA, CT and p are outside the "funnel" 
%                          =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  Trevor McDougall & Claire Roberts-Thomson.       [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (28th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  McDougall, T. J. and A. Klocker, 2010: An approximate geostrophic
%   streamfunction for use in density surfaces. Ocean Modelling, 32,
%   105-117.  
%    See section 8 of this paper.  
% 
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_geo_strf_dyn_height_pc:  Requires three inputs')
end %if
if ~(nargout == 2 | nargout == 3)
   error('gsw_geo_strf_dyn_height_pc:  Requires two or three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mdp,ndp] = size(delta_p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_dyn_height_pc: SA & CT need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_geo_strf_dyn_height_pc: There must be at least 2 values')
end

if (mdp == 1) & (ndp == 1)              % delta_p scalar - fill to size of SA
    error('gsw_geo_strf_dyn_height_pc: need more than one pressure')
elseif (ns == ndp) & (mdp == 1)         % delta_p is row vector,
    delta_p = delta_p(ones(1,ms), :);   % copy down each column.
elseif (ms == mdp) & (ndp == 1)          % delta_p is column vector,
    delta_p = delta_p(:,ones(1,ns));    % copy across each row.
elseif (ms == mdp) & (ns == ndp)
    % ok
else
    error('gsw_geo_strf_dyn_height_pc: Inputs array dimensions arguments do not agree')
end %if

transposed = 0;
if ms == 1  % row vector
   delta_p  =  delta_p(:);
   CT  =  CT(:);
   SA  =  SA(:);
   transposed = 1;
end %if

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

if ~isempty(find(delta_p < 0))
    error('gsw_geo_strf_dyn_height_pc: pressure must be monotonic')
end

p_deep = cumsum(delta_p); 
p_shallow = p_deep - delta_p;

in_funnel_shallow = gsw_infunnel(SA,CT,p_shallow);
in_funnel_deep = gsw_infunnel(SA,CT,p_deep);
in_funnel = in_funnel_shallow.*in_funnel_deep;

delta_h = gsw_enthalpy_diff_CT25(SA,CT,p_shallow,p_deep);

dyn_height_deep = -1*cumsum(delta_h);
%            This is Phi minus Phi_0 of Eqn. (3.32.2) of IOC et al. (2010).

p_mid = 0.5*(p_shallow  + p_deep);
delta_h_half = gsw_enthalpy_diff_CT25(SA,CT,p_mid,p_deep);

geo_strf_dyn_height_pc = gsw_enthalpy_SSO_0_CT25(p_mid) + ...
                           dyn_height_deep + delta_h_half;

%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly piecewise constant using 
% the computationally-efficient 25-term expression for density in terms of
% SA, CT and p. If one wanted to compute dynamic height anomaly with the 
% full TEOS-10 Gibbs function expression for density, the following lines 
% of code will enable this.
%
% delta_h = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep);
% dyn_height_deep = -1*cumsum(delta_h);
% p_mid = 0.5*(p_shallow  + p_deep);
% delta_h_half = gsw_enthalpy_diff_CT(SA,CT,p_mid,p_deep);
% SA_SO = 35.16504*ones(size(SA));
% CT_0 = zeros(size(CT));
% geo_strf_dyn_height_pc = gsw_enthalpy_CT(SA_SO,CT_0,p_mid) + ...
%                            dyn_height_deep + delta_h_half;
%
%---------------This is the end of the alternative code--------------------

if transposed
    geo_strf_dyn_height_pc = geo_strf_dyn_height_pc';
    p_mid = p_mid';
    in_funnel = in_funnel';
end

end
