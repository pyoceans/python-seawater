function [geo_strf_Cunningham, in_funnel] = gsw_geo_strf_Cunningham(SA,CT,p,interp_style)

% gsw_geo_strf_Cunningham             Cunningham geostrophic streamfunction
%==========================================================================
%
% USAGE:  
%  [geo_strf_Cunningham, in_funnel] = gsw_geo_strf_Cunningham(SA,CT,p,interp_style)
%
% DESCRIPTION:
%  Calculates the Cunningham geostrophic streamfunction (see Eqn. (3.29.2) 
%  of IOC et al. (2010)).  This is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned,
%  p, and the horizontal velocity at the sea surface.  This function 
%  calculates specific volume anomaly using the computationally efficient 
%  25-term expression for specific volume of McDougall et al. (2010).  
%  Under the default setting, this function evaluates the pressure integral
%  of specific volume using SA and CT “interploted” with respect to pressure
%  using a scheme based on the method of Reiniger and Ross (1968).  Our 
%  method uses a weighted mean of (i) values obtained from linear 
%  interpolation of the two nearest data points, and (ii) a linear 
%  extrapolation of the pairs of data above and below. This "curve fitting"
%  method resembles the use of cubic splines.  If the option “linear” is 
%  chosen, the function interpolates Absolute Salinity and Conservative 
%  Temperature linearly with presure in the vertical between “bottles”.
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  CT   =  Conservative Temperature                               [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%         ( ie. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  interp_style = interpolation technique.
%               = if nothing is entered the programme defaults to "curved" 
%                 interpolation between bottles in the vertical.
%               = if "linear" or "lin" is entered then the programme 
%                 interpolates linearly between bottles in the
%                 vertical.
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  geo_strf_Cunningham = Cunningham geostrophic streamfunction  [ m^2/s^2 ]
%
%  in_funnel           = 0, if SA, CT and p are outside the "funnel" 
%                      = 1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR:  
%   Trevor McDougall and Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  Cunningham, S.A., 2000: Circulation and volume flux of the North 
%   Atlantic using syoptic hydrographic data in a Bernoulli inverse.
%   J. Marine Res., 58, 1-35. 
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.29 of this TEOS-10 Manual. 
%
%  McDougall T.J., D.R. Jackett, P.M. Barker, C. Roberts-Thomson, R.
%   Feistel and R.W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3| nargin == 4)
   error('gsw_geo_strf_Cunningham:  Requires three or four inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_Cunningham: SA & CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    error('gsw_geo_strf_Cunningham : need more than one pressure'); 
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_geo_strf_Cunningham : Inputs array dimensions arguments do not agree')
end %if

if ~exist('interp_style','var')
    interp_style = 'curve';
elseif strcmpi('interp_style','linear') == 1 | strcmpi('interp_style','lin') == 1 |...
        strcmpi('interp_style','linaer') == 1 | strcmpi('interp_style','lnear') == 1
    interp_style = 'linear';
end

if ms == 1
    SA = SA';
    CT = CT';
    p = p';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

cp0 = 3991.86795711963;           % from Eqn. (3.3.3) of IOC et al. (2010).

delta_p = 1;  % This is the maximum distance between bottles in the vertical.

[dyn_height, in_funnel] = gsw_geo_strf_dyn_height(SA,CT,p,delta_p,interp_style);

geo_strf_Cunningham = dyn_height - gsw_enthalpy_SSO_0_CT25(p) + ...
                            gsw_enthalpy_CT25(SA,CT,p) - cp0*CT; 
                        
%--------------------------------------------------------------------------
% This function calculates the Cunningham streamfunction using the 
% computationally efficient 25-term expression for density in terms of SA, 
% CT and p.  If one wanted to compute this with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable 
% this.  Note that dynamic height will also need to be evaluated using the
% full Gibbs function.
%
%   SA_SO = 35.16504*ones(size(SA));
%   CT_0 = zeros(size(CT));
%   geo_strf_Cunningham = dyn_height - gsw_enthalpy_CT(SA_SO,CT_0,p) + ...
%                            gsw_enthalpy_CT(SA,CT,p) - cp0*CT; 
%
%---------------This is the end of the alternative code--------------------

if transposed
   geo_strf_Cunningham = geo_strf_Cunningham';
   in_funnel = in_funnel';
end %if

end
