function [pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA,CT)

% gsw_pt_second_derivatives     second derivatives of potential temperature
% =========================================================================
%
% USAGE:
%  [pt_SA_SA, pt_SA_CT, pt_CT_CT] = gsw_pt_second_derivatives(SA,CT)
%
% DESCRIPTION:
%  Calculates the following three second-order derivatives of potential
%  temperature (the regular potential temperature which has a reference
%  sea pressure of 0 dbar),
%   (1) pt_SA_SA, the second derivative with respect to Absolute Salinity
%       at constant Conservative Temperature,
%   (2) pt_SA_CT, the derivative with respect to Conservative Temperature
%       and Absolute Salinity, and
%   (3) pt_CT_CT, the second derivative with respect to Conservative
%       Temperature at constant Absolute Salinity.
%
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  CT  =   Conservative Temperature                               [ deg C ]
%
%  SA & CT need to have the same dimensions.
%
% OUTPUT:
%  pt_SA_SA  =  The second derivative of potential temperature (the
%               regular potential temperature which has reference sea
%               pressure of 0 dbar) with respect to Absolute Salinity
%               at constant Conservative Temperature.
%               pt_SA_SA has units of:                     [ K/((g/kg)^2) ]
%  pt_SA_CT  =  The derivative of potential temperature with respect
%               to Absolute Salinity and Conservative Temperature.
%               pt_SA_CT has units of:                        [ (g/kg)^-1 ]
%  pt_CT_CT  =  The second derivative of potential temperature (the
%               regular one with pr = 0 dbar) with respect to
%               Conservative Temperature at constant SA.
%               pt_CT_CT has units of:                             [ K^-1 ]
%
% AUTHOR:
%  Trevor McDougall and Paul Barker.        [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See Eqns. (A.12.9) and (A.12.10) of this TEOS-10 Manual.
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

if ~(nargin == 2)
   error('gsw_pt_second_derivatives:  Requires two inputs')
end %if

if ~(nargout == 3)
   error('gsw_pt_second_derivatives:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);

if (mt == 1) & (nt == 1)              % CT scalar - fill to size of SA
    CT = CT*ones(size(SA));
elseif (ns == nt) & (mt == 1)         % CT is row vector,
    CT = CT(ones(1,ms), :);              % copy down each column.
elseif (ms == mt) & (nt == 1)         % CT is column vector,
    CT = CT(:,ones(1,ns));               % copy across each row.
elseif (ms == mt) & (ns == nt)
    % ok
else
    error('gsw_pt_second_derivatives: The dimensions of CT do not agree')
end %if

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

dSA = 1e-3;                 % increment of Absolute Salinity is 0.001 g/kg.
SA_l = nan(size(SA));
inds_l = find(SA>=dSA);
if ~isempty(inds_l)
    SA_l(inds_l) = SA(inds_l) - dSA;
end
inds_l = find(SA<dSA);
if ~isempty(inds_l)
    SA_l(inds_l) = 0;
end
SA_u = SA + dSA;

[pt_SA_l, pt_CT_l] = gsw_pt_first_derivatives(SA_l,CT);
[pt_SA_u, pt_CT_u] = gsw_pt_first_derivatives(SA_u,CT);

pt_SA_SA = (pt_SA_u - pt_SA_l)./(SA_u - SA_l);
% pt_SA_CT = (pt_CT_u - pt_CT_l)./(SA_u - SA_l); % can calculate this either way;

dCT  = 1e-2;     % increment of Conservative Temperature is 0.01 degrees C;
CT_l = CT - dCT;
CT_u = CT + dCT;

[pt_SA_l, pt_CT_l] = gsw_pt_first_derivatives(SA,CT_l);
[pt_SA_u, pt_CT_u] = gsw_pt_first_derivatives(SA,CT_u);

pt_SA_CT = (pt_SA_u - pt_SA_l)./(CT_u - CT_l);
pt_CT_CT = (pt_CT_u - pt_CT_l)./(CT_u - CT_l);

if transposed
    pt_SA_CT = pt_SA_CT';
    pt_CT_CT = pt_CT_CT';
end

end
