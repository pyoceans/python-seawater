function [ntp_pt_vs_CT_ratio_CT25, in_funnel] = gsw_ntp_pt_vs_CT_ratio_CT25(SA,CT,p)

% gsw_ntp_pt_vs_CT_ratio_CT25               ratio of gradients of potential
%                             temperature and Conservative Temperature in a
%                           neutral tangent  plane (in a locally-referenced
%                              potential density surface)(25-term equation)
% =========================================================================
%
% USAGE:
%  [ntp_pt_vs_CT_ratio_CT25, in_funnel] = gsw_ntp_pt_vs_CT_ratio_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the ratio of the two-dimensional gradient of potential 
%  temperature versus that of Conservative Temperature, CT, along the  
%  neutral tangent plane. The potential temperature is the regular one  
%  which has a reference sea pressure of 0 dbar. Part of the calculation  
%  uses the computationally-efficient 25-term expression for density in   
%  terms of SA, CT and p (McDougall et al., 2010).  
%
% INPUT:
%  SA  =   Absolute Salinity                                       [ g/kg ]
%  CT  =   Conservative Temperature                               [ deg C ]
%  p   =   sea pressure                                            [ dbar ]
%          (ie. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  ntp_pt_vs_CT_ratio_CT25  =  The ratio of the spatial gradient of 
%                              potential temperature versus that of 
%                              Conservative Temperature in the 
%                              neutral tangent plane (ntp).    [ unitless ]
%  in_funnel  =  0, if SA, CT and p are outside the "funnel" 
%             =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker    [ gsw_help@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqn. (A.14.5) of this TEOS-10 Manual.   
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

if ~(nargin == 3)
   error('gsw_ntp_pt_vs_CT_ratio_CT25:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_ntp_pt_vs_CT_ratio_CT25: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_ntp_pt_vs_CT_ratio_CT25: Inputs array dimensions arguments do not agree')
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

in_funnel = gsw_infunnel(SA,CT,p);

[dummy, alpha_CT, beta_CT, dummy2] = gsw_rho_alpha_beta_CT25(SA,CT,p);

%--------------------------------------------------------------------------
% This function calculates the ntp_pt_vs_CT_ratio_CT25 using the 
% computationally-efficient 25-term expression for density in terms of SA, 
% CT and p. If one wanted to compute this with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable
% this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t = gsw_pt_from_t(SA,pt,pr0,p);
%    beta_CT = gsw_beta_const_CT(SA,t,p);
%    alpha_CT = gsw_alpha_wrt_CT(SA,t,p);
%
%--------- This is the end of the alternative code-------------------------

[pt_SA, pt_CT] = gsw_pt_first_derivatives(SA,CT);

ntp_pt_vs_CT_ratio_CT25 = pt_CT + pt_SA.*(alpha_CT./beta_CT);

if transposed
    ntp_pt_vs_CT_ratio_CT25 = ntp_pt_vs_CT_ratio_CT25';
    in_funnel = in_funnel';
end

end
