function [CT_SA, CT_pt] = gsw_CT_first_derivatives(SA,pt)

% gsw_CT_first_derivatives    first derivatives of Conservative Temperature
%==========================================================================
%
% USAGE:
%  [CT_SA, CT_pt] = gsw_CT_first_derivatives(SA,pt)
%
% DESCRIPTION:
%  Calculates the following two derivatives of Conservative Temperature
%  (1) CT_SA, the derivative with respect to Absolute Salinity at
%      constant potential temperature (with pr = 0 dbar), and
%   2) CT_pt, the derivative with respect to potential temperature
%      (the regular potential temperature which is referenced to 0 dbar)
%      at constant Absolute Salinity.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  pt  =  potential temperature (ITS-90)                          [ deg C ]
%         (whose reference pressure is 0 dbar)
%
%  SA & pt need to have the same dimensions.
%
% OUTPUT:
%  CT_SA  =  The derivative of Conservative Temperature with respect to
%            Absolute Salinity at constant potential temperature
%            (the regular potential temperature which has reference
%            sea pressure of 0 dbar). The CT_SA output has units of:
%                                                               [ K/(g/kg)]
%  CT_pt  =  The derivative of Conservative Temperature with respect to
%            potential temperature (the regular one with pr = 0 dbar)
%            at constant SA. CT_pt is dimensionless.           [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker   [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.3) and (A.12.9a,b) of this TEOS-10 Manual.
%
% McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%  Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
%  expression for the density of seawater in terms of Conservative
%  Temperature, and related properties of seawater.  To be submitted
%  to Ocean Science Discussions.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_CT_first_derivatives:  Requires two inputs')
end %if

if ~(nargout == 2)
   error('gsw_CT_first_derivatives:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(pt);

if (mt ~= ms | nt ~= ns)
    error('gsw_CT_first_derivatives: SA and t must have same dimensions')
end

if ms == 1
    SA = SA';
    pt = pt';
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
abs_pt = 273.15 + pt;

CT_SA = (gsw_gibbs(n1,n0,n0,SA,pt,pr0) -...
               abs_pt.*gsw_gibbs(n1,n1,n0,SA,pt,pr0))./cp0;

CT_pt = - (abs_pt.*gsw_gibbs(n0,n2,n0,SA,pt,pr0))./cp0;

if transposed
    CT_SA = CT_SA';
    CT_pt = CT_pt';
end

end
