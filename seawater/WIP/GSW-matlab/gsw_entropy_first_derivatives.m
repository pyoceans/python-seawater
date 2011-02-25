function [eta_SA, eta_CT] = gsw_entropy_first_derivatives(SA,CT)

% gsw_entropy_first_derivatives                first derivatives of entropy
% =========================================================================
%
% USAGE:
%  [eta_SA, eta_CT] = gsw_eta_first_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following two partial derivatives of specific entropy
%  (eta)
%   (1) eta_SA, the derivative with respect to Absolute Salinity at
%       constant Conservative Temperature, and
%   (2) eta_CT, the derivative with respect to Conservative Temperature at
%       constant Absolute Salinity.
%
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  CT  =   Conservative Temperature                               [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  eta_SA =  The derivative of specific entropy with respect to
%            Absolute Salinity (in units of g/kg) at constant Conservative
%            Temperature.  The eta_SA output has units of:
%                                         [ J/(kg K(g/kg))]  or [ J/(g K) ]
%  eta_CT =  The derivative of specific entropy with respect to
%            Conservative Temperature at constant Absolute Salinity.
%            The eta_CT output has units of:                 [ J/(kg K^2) ]
%
% AUTHOR:
%  Trevor McDougall.  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (21th August, 2010).
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.8) and (P.14a,c) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_entropy_first_derivatives:  Requires two inputs')
end %if

if ~(nargout == 2)
   error('gsw_entropy_first_derivatives:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_entropy_first_derivatives: SA and CT must have same dimensions')
end
if ms == 1
    SA = SA';
    CT = CT';
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
pt = gsw_pt_from_CT(SA,CT);

eta_SA = -(gsw_gibbs(n1,n0,n0,SA,pt,pr0))./(273.15 + pt);

eta_CT = (cp0*ones(size(pt)))./(273.15 + pt);

if transposed
    eta_SA = eta_SA';
    eta_CT = eta_CT';
end

end
