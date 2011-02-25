function [h_SA, h_CT, h_P] = gsw_enthalpy_first_derivatives(SA,CT,p)

% gsw_enthalpy_first_derivatives              first derivatives of enthalpy
%==========================================================================
%
% USAGE:
%  [h_SA, h_CT, h_P] = gsw_enthalpy_first_derivatives(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the following three derivatives of specific enthalpy (h)
%   (1) h_SA, the derivative with respect to Absolute Salinity at
%       constant CT and p, and
%   (2) h_CT, derivative with respect to CT at constant SA and p.
%   (3) h_P, derivative with respect to pressure (in Pa) at constant
%       SA and CT.
%
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  CT  =   Conservative Temperature                               [ deg C ]
%  p   =   sea pressure                                            [ dbar ]
%          (i.e. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  h_SA   =  The first derivative of specific enthalpy with respect to
%            Absolute Salinity at constant CT and p.
%                                            [ J/(kg (g/kg))]  i.e. [ J/g ]
%  h_CT   =  The first derivative of specific enthalpy with respect to
%            CT at constant SA and p.                          [ J/(kg K) ]
%  h_P    =  The first partial derivative of specific enthalpy with
%            respect to pressure (in Pa) at fixed SA and CT.  Note that
%            h_P is specific volume (1/rho).
%
% AUTHOR:
%  Trevor McDougall.  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (17th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.11.18), (A.11.15) and (A.11.12) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_enthalpy_first_derivatives:  Requires three inputs')
end %if

if ~(nargout == 3)
   error('gsw_enthalpy_first_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms ~= mt | ns ~= nt )
    error('gsw_enthalpy_first_derivatives: SA and CT do not have the same dimensions')
end %if

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA.
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_enthalpy_first_derivatives: The dimensions of p do not agree')
end %if

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
n0 = 0;
n1 = 1;
pr0 = zeros(size(SA));

pt0 = gsw_pt_from_CT(SA,CT);
t   = gsw_pt_from_t(SA,pt0,pr0,p);
temp_ratio = (273.15 + t)./(273.15 + pt0);

h_CT = cp0.*temp_ratio;         % from Eqn. (A.11.15) of IOC et al. (2010).

h_SA = gsw_gibbs(n1,n0,n0,SA,t,p) - ...
            temp_ratio.*gsw_gibbs(n1,n0,n0,SA,pt0,pr0);
                                % from Eqn. (A.11.18) of IOC et al. (2010).
h_P = gsw_gibbs(n0,n0,n1,SA,t,p);
                                % from Eqn. (A.11.12) of IOC et al. (2010).
if transposed
    h_CT = h_CT';
    h_SA = h_SA';
    h_P = h_P';
end

end