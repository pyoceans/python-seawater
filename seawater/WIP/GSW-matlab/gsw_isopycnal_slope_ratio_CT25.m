function [isopycnal_slope_ratio_CT25, in_funnel] = gsw_isopycnal_slope_ratio_CT25(SA,CT,p,pr)

% gsw_isopycnal_slope_ratio_CT25          ratio of the slopes of isopycnals
%                                         on the SA-CT diagram for p and pr
%                                                        (25-term equation)
% =========================================================================
%
% USAGE:
% [isopycnal_slope_ratio_CT25, in_funnel] = gsw_isopycnal_slope_ratio_CT25(SA,CT,p,pr)
%
% DESCRIPTION:
%  Calculates the ratio of alpha_CT/beta_CT at pressure p to that at 
%  pressure pr.  This function uses the computationally-efficient 25-term 
%  expression for density in terms of SA, CT & p (McDougall et al., 2010).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%  pr  =  reference pressure                                       [ dbar ]
%         (ie. absolute reference pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p & pr may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN
%
% OUTPUT:
%  isopycnal_slope_ratio_CT25  
%                         =  The ratio of alpha_CT/beta_CT evaluated at 
%                            pressure p to that at pressure pr.  
%                                                              [ unitless ]
%  in_funnel              =  0, if SA, CT and p are outside the "funnel" 
%                         =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  Trevor McDougall, Paul Barker & David Jackett [ gsw_help@csiro.au ]
%      
% VERSION NUMBER: 2.0 (26th August, 2010).
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.  
%    See Eqn. (3.17.2) of this TEOS-10 Manual.   
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

if ~(nargin == 4)
   error('gsw_isopycnal_slope_ratio_CT25:  Requires four inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[mpr,npr] = size(pr);

if (mt ~= ms | nt ~= ns)
    error('gsw_isopycnal_slope_ratio_CT25: SA and CT must have same dimensions')
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
    error('gsw_isopycnal_slope_ratio_CT25: Inputs array dimensions arguments do not agree')
end %if

if (mpr == 1) & (npr == 1)              % pr scalar - fill to size of SA
    pr = pr*ones(size(SA));
elseif (ns == npr) & (mpr == 1)         % pr is row vector,
    pr = pr(ones(1,ms), :);              % copy down each column.
elseif (ms == mpr) & (npr == 1)         % pr is column vector,
    pr = pr(:,ones(1,ns));               % copy across each row.
elseif (ms == mpr) & (ns == npr)
    % ok
else
    error('gsw_isopycnal_slope_ratio_CT25: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    CT = CT';
    p = p';
    pr = pr';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = gsw_infunnel(SA,CT,p);

[dummy, alpha_CT, beta_CT, dummy2] = gsw_rho_alpha_beta_CT25(SA,CT,p);
[dummy, alpha_CT_pr, beta_CT_pr, dummy2] = gsw_rho_alpha_beta_CT25(SA,CT,pr);

%-------------------------------------------------------------------------
% This function calculates isopycnal_slope_ratio_CT25 using the 
% computationally-efficient 25-term expression for density as a function of
% SA, CT and p. If one wanted to compute this with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable 
% this.
%  
%     pt = gsw_pt_from_CT(SA,CT);
%     pr0 = zeros(size(SA)); 
%     t = gsw_pt_from_t(SA,pt,pr0,p);
%     beta_CT = gsw_beta_const_CT(SA,t,p);
%     alpha_CT = gsw_alpha_wrt_CT(SA,t,p);
%     tr = gsw_pt_from_t(SA,pt,pr0,pr);
%     beta_CT_pr = gsw_beta_const_CT(SA,tr,pr);
%     alpha_CT_pr = gsw_alpha_wrt_CT(SA,tr,pr);
%
%--------------This is the end of the alternative code---------------------

isopycnal_slope_ratio_CT25 = nan(size(SA));
[I] = find(alpha_CT_pr ~= 0);
if ~isempty(I)
    isopycnal_slope_ratio_CT25 (I) = (alpha_CT(I).*beta_CT_pr(I))./ ...
        (alpha_CT_pr(I).*beta_CT(I));
end

if transposed
    isopycnal_slope_ratio_CT25 = isopycnal_slope_ratio_CT25';
    in_funnel = in_funnel';
end

end
