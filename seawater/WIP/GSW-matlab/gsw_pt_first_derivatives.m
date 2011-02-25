function [pt_SA, pt_CT] = gsw_pt_first_derivatives(SA,CT)

% gsw_pt_first_derivatives       first derivatives of potential temperature
% =========================================================================
%
% USAGE:
%   [pt_SA, pt_CT] = gsw_pt_first_derivatives(SA,CT)
%
% DESCRIPTION:
% Calculates the following two partial derivatives of potential temperature
% (the regular potential temperature whose reference sea pressure is 0 dbar)
%   (1) pt_SA, the derivative with respect to Absolute Salinity at
%       constant Conservative Temperature, and
%   (2) pt_CT, the derivative with respect to Conservative Temperature at
%       constant Absolute Salinity.
%
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  CT  =   Conservative Temperature                               [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt_SA =  The derivative of potential temperature with respect to
%           Absolute Salinity at constant Conservative Temperature.
%                                                               [ K/(g/kg)]
%  pt_CT =  The derivative of potential temperature with respect to
%           Conservative Temperature at constant Absolute Salinity.
%           pt_CT is dimensionless.                            [ unitless ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker [ help_gsw@csiro.au ]
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.6), (A.12.3), (P.6) and (P.8) of this TEOS-10 Manual.
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term
%   expression for the density of seawater in terms of Conservative
%   Temperature, and related properties of seawater.  To be submitted
%   to Ocean Science Discussions.
%
%  This software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if (nargin == 3)
   error('gsw_pt_first_derivatives:  Requires three inputs')
end %if

if ~(nargout == 2)
   error('gsw_pt_first_derivatives:  Requires two outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt ~= ms | nt ~= ns)
    error('gsw_pt_first_derivatives: SA and CT must have same dimensions')
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
n2 = 2;
pr0 = zeros(size(SA));
pt = gsw_pt_from_CT(SA,CT);
abs_pt = (273.15 + pt);

CT_SA = (gsw_gibbs(n1,n0,n0,SA,pt,pr0) -...
                abs_pt.*gsw_gibbs(n1,n1,n0,SA,pt,pr0))./cp0;

CT_pt = - (abs_pt.*gsw_gibbs(n0,n2,n0,SA,pt,pr0))./cp0;

pt_SA = - CT_SA./CT_pt;

pt_CT = ones(size(CT_pt))./CT_pt;

if transposed
    pt_SA = pt_SA';
    pt_CT = pt_CT';
end

end
