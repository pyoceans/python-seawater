function [Tu, Rsubrho, p_mid, in_funnel] = gsw_Turner_Rsubrho_CT25(SA,CT,p)

% gsw_Turner_Rsubrho_CT25                            Turner angle & Rsubrho
%==========================================================================
% 
% USAGE:  
% [Tu, Rsubrho, p_mid, in_funnel] = gsw_Turner_Rsubrho_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the Turner angle and the Rsubrho as a function of pressure 
%  down a vertical water column.  These quantities express the relative 
%  contributions of the vertical gradients of Conservative Temperature 
%  and Absolute Salinity to the vertical stability (the square of the 
%  Brunt-Vaisala Frequency squared, N^2).  Tu and Rsubrho are evaluated at 
%  the mid pressure between the individual data points in the vertical.  
%  This function uses computationally-efficient 25-term expression for 
%  density in terms of SA, CT and p (McDougall et al., 2010). 
%
% INPUT:  
%   SA  =  Absolute Salinity                                       [ g/kg ]
%   CT  =  Conservative Temperature                               [ deg C ]
%   p   =  sea pressure                                            [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%   SA & CT need to have the same dimensions, 
%   p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%   Tu       =    Turner angle, on the same (M-1)xN grid as p_mid.
%                 Turner angle has units of:        [ degrees of rotation ]
%   Rsubrho  =    Stability Ratio, on the same (M-1)xN grid as p_mid.
%                 Rsubrho is dimensionless.                    [ unitless ]
%   p_mid    =    mid pressure between the indivual points of the p grid. 
%                 That is, p_mid is on a (M-1)xN grid in the vertical.  
%                 p_mid has units of:                              [ dbar ]
%   in_funnel  =  0, if (SA, CT and p) are outside the "funnel" 
%              =  1, if (SA, CT and p) are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR:  
%        Trevor McDougall & Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.15.1) and (3.16.1) of this TEOS-10 Manual. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_Turner_Rsubrho_CT25:  Requires three inputs')
end %if

if ~(nargout == 3 | nargout == 4)
   error('gsw_Turner_Rsubrho_CT25:  Requires three or four outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns )
    error('gsw_Turner_Rsubrho_CT25: SA and CT must have same dimensions')
end

if (ms*ns == 1)
    error('There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p is a scalar - must be two bottles
    error('There must be at least 2 pressure values')
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_Turner_Rsubrho_CT25: Inputs array dimensions arguments do not agree')
end %if

[mp,np] = size(p);

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

% These few lines ensure that SA is non-negative.
[I_neg_SA] = find(SA < 0);
if ~isempty(I_neg_SA)
    SA(I_neg_SA) = 0;
end
in_funnel = gsw_infunnel(SA,CT,p);

Ishallow = 1:(mp-1);
Ideep = 2:mp;
p_mid = (p(Ishallow,:) + p(Ideep,:))/2;
SA_mid = (SA(Ishallow,:) + SA(Ideep,:))/2;
CT_mid = (CT(Ishallow,:) + CT(Ideep,:))/2;

dSA = SA(Ishallow,:) - SA(Ideep,:);
dCT = CT(Ishallow,:) - CT(Ideep,:);

[dummy, alpha, beta, dummy2] = gsw_rho_alpha_beta_CT25(SA_mid,CT_mid,p_mid);

%--------------------------------------------------------------------------
% This function evaluates Tu and Rsubrho using the computationally-efficient
% 25-term expression for density in terms of SA, CT and p. If one wanted to
% compute Tu and Rsubrho using the full TEOS-10 Gibbs function expression 
% for density, the following lines of code would do that.  
%
%    pt_mid = gsw_pt_from_CT(SA_mid,CT_mid);
%    pr0 = zeros(size(SA_mid)); 
%    t_mid = gsw_pt_from_t(SA_mid,pt_mid,pr0,p_mid);
%    beta = gsw_beta_const_CT(SA_mid,t_mid,p_mid);
%    alpha = gsw_alpha_wrt_CT(SA_mid,t_mid,p_mid);
%
% --------------This is the end of the alternative code--------------------

Tu = atan2((alpha.*dCT + beta.*dSA),(alpha.*dCT - beta.*dSA));
Tu = Tu.*(180/pi);

Rsubrho = nan(size(dSA));
[Inz] = find(dSA ~= 0);
if ~isempty(Inz)
    Rsubrho(Inz) = (alpha(Inz).*dCT(Inz))./(beta(Inz).*dSA(Inz));
end

if transposed
    Tu      = Tu';
    Rsubrho = Rsubrho';
    p_mid   = p_mid';
    in_funnel = in_funnel';
end

end
