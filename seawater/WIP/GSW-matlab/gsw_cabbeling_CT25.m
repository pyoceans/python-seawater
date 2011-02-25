function [cabbeling_CT25, in_funnel] = gsw_cabbeling_CT25(SA,CT,p)

% gsw_cabbeling_CT25                                  cabbeling coefficient 
%                                                        (25-term equation)
%==========================================================================
%
% USAGE:  
%  [cabbeling_CT25, in_funnel] = gsw_cabbeling_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the cabbeling coefficient of seawater with respect to  
%  Conservative Temperature. This function uses the computationally-
%  efficient 25-term expression for density in terms of SA, CT and p
%  (McDougall et al., 2010)
%   
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (ie. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  cabbeling_CT25  =  Cabbeling coefficient with respect to      
%                     Conservative Temperature.                 [ 1/(K^2) ]
%  in_funnel       =  0, if SA, CT and p are outside the "funnel" 
%                  =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%   David Jackett, Trevor McDougall and Paul Barker  [ help_gsw@csiro.au ]   
%
% VERSION NUMBER: 2.0 (24th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.9.2) and (P.4) of this TEOS-10 manual.
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, 
%   R. Feistel and R. W. Hallberg, 2010:  A computationally efficient 
%   25-term expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted to 
%   Ocean Science Discussions.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_cabbeling_CT25:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_cabbeling_CT25: SA and t must have same dimensions')
end

if (mp == 1) & (np == 1)              % p is a scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_cabbeling_CT25: Inputs array dimensions arguments do not agree')
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

dCT = 1e-3;          % increment in Conservative Temperature is 1e-3 deg C.
CT_l = CT - dCT;
CT_u = CT + dCT;

[dummy,alpha_CT,beta_CT,flag] = gsw_rho_alpha_beta_CT25(SA,CT,p);
[dummy,alpha_CT_u,dummy2,flag] = gsw_rho_alpha_beta_CT25(SA,CT_u,p);
[dummy,alpha_CT_l,dummy2,flag] = gsw_rho_alpha_beta_CT25(SA,CT_l,p);

ratio_CT = alpha_CT./beta_CT;
alpha_CT_CT = (alpha_CT_u - alpha_CT_l)./(CT_u-CT_l);

dSA = 1e-3;                  % increment in Absolute Salinity is 1e-3 g/kg.
inds_l = find(SA>=dSA);
SA_l = nan(size(SA));
if ~isempty(inds_l)
    SA_l(inds_l) = SA(inds_l) - dSA;
end
inds_l = find(SA<dSA);
if ~isempty(inds_l)
    SA_l(inds_l) = 0;
end
SA_u = SA + dSA;  

[dummy,alpha_CT_u,beta_CT_u,flag] = gsw_rho_alpha_beta_CT25(SA_u,CT,p);
[dummy,alpha_CT_l,beta_CT_l,flag] = gsw_rho_alpha_beta_CT25(SA_u,CT,p);

alpha_CT_SA = (alpha_CT_u - alpha_CT_l)./(SA_u-SA_l);
beta_CT_SA = (beta_CT_u - beta_CT_l)./(SA_u-SA_l);
cabbeling_CT25 = alpha_CT_CT + ratio_CT.*(2.*alpha_CT_SA - ratio_CT.*beta_CT_SA);

%--------------------------------------------------------------------------
% This function calculates cabbeling_CT25 using the computationally-efficient
% 25-term expression for density in terms of SA, CT and p. If one wanted to
% compute cabbeling_CT25 with the full TEOS-10 Gibbs function expression 
% for density, the following lines of code will enable this.
%
%   pr0 = zeros(size(p)); 
%   pt = gsw_pt_from_CT(SA,CT);    
%   dpt = 1e-3;          % increment in potential temperature is 1e-3 deg C
%   pt_l = pt - dpt; 
%   pt_u = pt + dpt; 
%   CT_l = gsw_CT_from_pt(SA,pt_l);
%   CT_u = gsw_CT_from_pt(SA,pt_u);  
%   t_l = gsw_pt_from_t(SA,pt_l,pr0,p); 
%   t_u = gsw_pt_from_t(SA,pt_u,pr0,p);
%   t = 0.5*(t_l + t_u);
%   alpha_CT = gsw_alpha_wrt_CT(SA,t,p); 
%   beta_CT = gsw_beta_const_CT(SA,t,p); 
%   ratio_CT = alpha_CT./beta_CT;
%   alpha_CT_CT = (gsw_alpha_wrt_CT(SA,t_u,p)-gsw_alpha_wrt_CT(SA,t_l,p))./(CT_u-CT_l);
%   dSA = 1e-3;                %increment in Absolute Salinity is 1e-3 g/kg
%   SA_l = nan(size(SA));
%   inds_l = find(SA>=dSA);
%   if ~isempty(inds_l)   
%     SA_l(inds_l) = SA(inds_l)-dSA;
%   end
%   inds_l = find(SA<dSA);
%   if ~isempty(inds_l)   
%     SA_l(inds_l) = 0; 
%   end
%   SA_u = SA+dSA;  
%   alpha_CT_SA  = (gsw_alpha_wrt_CT(SA_u,t,p)-gsw_alpha_wrt_CT(SA_l,t,p))./(SA_u-SA_l);
%   beta_CT_SA   = (gsw_beta_const_CT(SA_u,t,p)-gsw_beta_const_CT(SA_l,t,p))./(SA_u-SA_l);
%   cabbeling_CT = alpha_CT_CT + ratio_CT.*(2.*alpha_CT_SA - ratio_CT.*beta_CT_SA);
%
%---------------This is the end of the alternative code--------------------

if transposed
    cabbeling_CT25 = cabbeling_CT25';
    in_funnel = in_funnel';
end

end
