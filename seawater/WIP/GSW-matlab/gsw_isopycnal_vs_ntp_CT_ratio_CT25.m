function [G_CT, p_mid, in_funnel] = gsw_isopycnal_vs_ntp_CT_ratio_CT25(SA,CT,p,pr)

% gsw_isopycnal_vs_ntp_CT_ratio_CT25               ratio of the gradient of
%                   Conservative Temperature in a potential density surface
%                     to that in a neutral tangent plane (i.e. in a locally
%                  referenced potential density surface) (25-term equation)
%==========================================================================
% 
% USAGE:  
%  [G_CT, p_mid, in_funnel] = gsw_isopycnal_vs_ntp_CT_ratio_CT25(SA,CT,p,pr)
%
% DESCRIPTION:
%  Calculates the ratio of the two-dimensional gradient of Conservative
%  Temperature in a potential density surface (with reference sea pressure
%  (pr)) versus that in the neutral tangent plane (see Eqns. (3.17.3) and
%  (3.17.4) of IOC et al. (2010)).  This ratio has been called the
%  "isopycnal Conservative Temperature gradient ratio".  This ratio is
%  evaluated at the mid pressure between the individual data points in the
%  vertical. The reference sea pressure of the potential density surface
%  must have a constant value.  This function uses the computationally
%  efficient 25-term expression for density in terms of SA, CT and p
%  (McDougall et al., 2010).  
%
% INPUT:  
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%  pr  =  reference sea pressure of the potential density surface
%         (ie. absolute reference pressure - 10.1325 dbar)         [ dbar ]
%
%  SA & CT need to have the same dimensions.
%  p & pr may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT 
%  are MxN
%
% OUTPUT:
%  G_CT      =  the ratio of the gradient of CT in a potential density 
%               surface to that in a neutral tangent plane.  
%               G_CT is output on the same vertical (M-1)xN grid 
%               as p_mid, where M & N are the dimensions of SA.
%               G_CT is dimensionless.                         [ unitless ]
%  p_mid     =  mid pressure between the individual points of the p grid. 
%               That is, p_mid is on a (M-1)xN grid.  
%               p_mid has units of:                                [ dbar ]
%  in_funnel =  0, if SA, CT and p are outside the "funnel" 
%            =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR:  
%  Trevor McDougall, Paul Barker & David Jackett [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.17.3) and (3.17.4) of this TEOS-10 Manual. 
%
%  McDougall, T. J., 1987: Neutral surfaces. Journal of Physical 
%   Oceanography, 17, 1950-1964.  See Eqn. (29) of this paper.  
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

if ~(nargin == 3 | nargin == 4)
   error('gsw_isopycnal_vs_ntp_CT_ratio_CT25:  Requires three or four inputs')
end %if

if ~(nargout == 2 | nargout == 3)
   error('gsw_isopycnal_vs_ntp_CT_ratio_CT25:  Requires two or three outputs')
end

if nargin == 3
%    Assume reference pressure is 0 dbar.
  pr = 0;
end %if

if ~isscalar(unique(pr))
    error('gsw_IPV_vs_fNsquared_ratio_CT25: The reference pressures differ, they should be unique')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('SA and CT must have same dimensions, p can be a vector')
end

if (ms*ns == 1)
    error('There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p scalar - must be two bottles
    error('There must be at least 2 values')
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_isopycnal_vs_ntp_CT_ratio_CT25: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    CT = CT';
    p = p';
    [mp,np] = size(p);
    transposed = 1;
else
    transposed = 0;
end

pr = unique(pr)*ones(mp-1,np);               %resize the reference pressure 

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = gsw_infunnel(SA,CT,p);

Ishallow = 1:(mp-1);
Ideep = 2:mp;
p_mid = 0.5*(p(Ishallow,:) + p(Ideep,:));
SA_mid = 0.5*(SA(Ishallow,:) + SA(Ideep,:));
CT_mid = 0.5*(CT(Ishallow,:) + CT(Ideep,:));

dSA = SA(Ishallow,:) - SA(Ideep,:);
dCT = CT(Ishallow,:) - CT(Ideep,:);

[dummy, alpha, beta, dummy2] = gsw_rho_alpha_beta_CT25(SA_mid,CT_mid,p_mid);
[dummy, alpha_pr, beta_pr, dummy2] = gsw_rho_alpha_beta_CT25(SA_mid,CT_mid,pr);

%--------------------------------------------------------------------------
% This function calculates G_CT using the computationally-efficient 
% 25-term expression for density as a function of SA, CT and p. If one 
% wanted to compute this with the full TEOS-10 Gibbs function expression 
% for density, the following lines of code will enable this.
%
%    pt_mid = gsw_pt_from_CT(SA_mid,CT_mid);
%    pr0 = zeros(size(SA_mid)); 
%    t_mid = gsw_pt_from_t(SA_mid,pt_mid,pr0,p_mid);
%    beta = gsw_beta_const_CT(SA_mid,t_mid,p_mid);
%    alpha = gsw_alpha_wrt_CT(SA_mid,t_mid,p_mid);
%    beta_pr  = gsw_beta_const_CT(SA_mid,t_mid,pr);
%    alpha_pr = gsw_alpha_wrt_CT(SA_mid,t_mid,pr);
%
%------------This is the end of the alternative code----------------------- 

anum   = dCT.*alpha./beta - dSA;
adenom = dCT.*alpha_pr./beta_pr - dSA;

G_CT = nan(size(SA_mid));
[I] = find(adenom ~= 0);
if ~isempty(I)
    G_CT(I) = anum(I)./adenom(I);
end

if transposed
    G_CT    = G_CT';
    p_mid   = p_mid';
    in_funnel = in_funnel';
end

end
