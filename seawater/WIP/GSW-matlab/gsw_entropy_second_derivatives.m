function [eta_SA_SA, eta_SA_CT, eta_CT_CT] = gsw_entropy_second_derivatives(SA,CT)

% gsw_entropy_second_derivatives              second derivatives of entropy
% =========================================================================
%
% USAGE:
%  [eta_SA_SA, eta_SA_CT, eta_CT_CT] = gsw_entropy_second_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following three second-order partial derivatives of
%  specific entropy (eta)
%   (1) eta_SA_SA, the second derivative with respect to Absolute
%       Salinity at constant Conservative Temperature, and
%   (2) eta_SA_CT, the derivative with respect to Absolute Salinity and
%       Conservative Temperature.
%   (3) eta_CT_CT, the second derivative with respect to Conservative
%       Temperature at constant Absolute Salinity.
%
% INPUT:
%  SA   =   Absolute Salinity                                      [ g/kg ]
%  CT   =   Conservative Temperature                              [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  eta_SA_SA =  The second derivative of specific entropy with respect
%               to Absolute Salinity (in units of g/kg) at constant
%               Conservative Temperature.  eta_SA_SA has units of:
%                                                       [ J/(kg K(g/kg)^2)]
%  eta_SA_CT =  The second derivative of specific entropy with respect
%               to Conservative Temperature at constant Absolute
%               Salinity. eta_SA_CT has units of:     [ J/(kg (g/kg) K^2) ]
%
%  eta_CT_CT =  The second derivative of specific entropy with respect
%               to Conservative Temperature at constant Absolute
%               Salinity.  eta_CT_CT has units of:           [ J/(kg K^3) ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (P.14b) and (P.15a,b) of this TEOS-10 Manual.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_entropy_second_derivatives:  Requires two inputs')
end %if

if ~(nargout == 3)
   error('gsw_entropy_second_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_entropy_second_derivatives: SA and CT must have same dimensions')
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

cp0 = 3991.86795711963;             % from Eqn. 3.3.3 of IOC et al. (2010).
n0 = 0;
n1 = 1;
n2 = 2;
pr0 = zeros(size(SA));
pt = gsw_pt_from_CT(SA,CT);
abs_pt = 273.15 + pt;

CT_SA = (gsw_gibbs(n1,n0,n0,SA,pt,pr0) - ...
               (abs_pt.*gsw_gibbs(n1,n1,n0,SA,pt,pr0)))./cp0;

CT_pt = -(abs_pt.*gsw_gibbs(n0,n2,n0,SA,pt,pr0))./cp0;

eta_CT_CT = -(cp0*ones(size(SA)))./(CT_pt.*abs_pt.*abs_pt);

eta_SA_CT = - CT_SA.*eta_CT_CT;

eta_SA_SA = - gsw_gibbs(n2,n0,n0,SA,pt,pr0)./abs_pt - CT_SA.*eta_SA_CT;

if transposed
    eta_CT_CT = eta_CT_CT';
    eta_SA_CT = eta_SA_CT';
    eta_SA_SA = eta_SA_SA';
end

end
