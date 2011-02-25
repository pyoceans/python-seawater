function [thermobaric_CT25, in_funnel] = gsw_thermobaric_CT25(SA,CT,p)

% gsw_thermobaric_CT25           thermobaric coefficient (25-term equation)
%==========================================================================
%
% USAGE:  
%  [thermobaric_CT25, in_funnel] = gsw_thermobaric_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the thermobaric coefficient of seawater with respect to
%  Conservative Temperature.  This routine calculates rho from the 
%  computationally-efficient 25-term expression for density in terms of
%  SA, CT and p (McDougall et al., 2010).
%  
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( ie. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  thermobaric_CT25  =  Thermobaric coefficient with           [ 1/(K Pa) ] 
%                       respect to Conservative Temperature.           
%  Note. The pressure derivative is taken with respect to
%    pressure in Pa not dbar.
%
%  in_funnel         =  0, if SA, CT and p are outside the "funnel" 
%                    =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (24th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.8.2) and (P.2) of this TEOS-10 manual.
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
   error('gsw_thermobaric_CT25:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_thermobaric_CT25: SA and CT must have same dimensions')
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
    error('gsw_thermobaric_CT25: Inputs array dimensions arguments do not agree')
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

db2Pa = 1e4;
dp = 1e-1;                        % pressure increment is 1e-1 dbar (10 Pa)

if p==0, p = p+dp; end
inds = find(p>=dp); 
p_u = zeros(size(p));
p_l = dp*ones(size(p));
if ~isempty(inds)
    p_u(inds) = p(inds)-dp;
    p_l(inds) = p(inds)+dp;
end

[dummy,alpha_CT,beta_CT,flag] = gsw_rho_alpha_beta_CT25(SA,CT,p);
[dummy,alpha_CT_u,beta_CT_u,flag] = gsw_rho_alpha_beta_CT25(SA,CT,p_u);
[dummy,alpha_CT_l,beta_CT_l,flag] = gsw_rho_alpha_beta_CT25(SA,CT,p_l);

alpha_CT_p = (alpha_CT_u - alpha_CT_l)./(p_u-p_l);
beta_CT_p  = (beta_CT_u - beta_CT_l)./(p_u-p_l);

thermobaric_CT = alpha_CT_p - (alpha_CT./beta_CT).*beta_CT_p;
thermobaric_CT25 = thermobaric_CT./db2Pa;         % To have units of 1/(K Pa)

%--------------------------------------------------------------------------
% This function calculates thermobaric_CT using the computationally-efficient 
% 25-term expression for density in terms of SA, CT and p. If one wanted to
% compute thermobaric_CT with the full TEOS-10 Gibbs function expression 
% for density, the following lines of code will do this.
%
%  pr0 = zeros(size(p)); 
%  pt = gsw_pt_from_CT(SA,CT);
%  t_l = gsw_pt_from_t(SA,pt,pr0,p_l);   
%  t_u = gsw_pt_from_t(SA,pt,pr0,p_u);
%  t = 0.5*(t_l + t_u);
%  alpha_CT = gsw_alpha_wrt_CT(SA,t,p);
%  beta_CT = gsw_beta_const_CT(SA,t,p);
%  alpha_CT_p = (gsw_alpha_wrt_CT(SA,t_u,p_u)-gsw_alpha_wrt_CT(SA,t_l,p_l))./(p_u-p_l);
%  beta_CT_p = (gsw_beta_const_CT(SA,t_u,p_u)-gsw_beta_const_CT(SA,t_l,p_l))./(p_u-p_l);
%  thermobaric_CT = alpha_CT_p - (alpha_CT./beta_CT).*beta_CT_p;
%  thermobaric_CT = thermobaric_CT./db2Pa;      % To have units of 1/(K Pa)
%
%----------------This is the end of the alternative code-------------------

if transposed
    thermobaric_CT25 = thermobaric_CT25';
    in_funnel = in_funnel';
end

end
