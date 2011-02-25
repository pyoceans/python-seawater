function [rho_CT25, in_funnel] = gsw_rho_CT25(SA,CT,p)

% gsw_rho_CT25                           In-situ density (25-term equation)
%==========================================================================
% 
% USAGE:  
%  rho_CT25 = gsw_rho_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density from Absolute Salinity and Conservative 
%  Temperature, using the computationally-efficient 25-term expression 
%  for density in terms of SA, CT and p (McDougall et al., 2010).
%
%  Note that potential density with respect to reference pressure, pr, is
%  obtained by calling this function with the pressure argument being pr
%  (i.e. "gsw_rho_CT25(SA,CT,pr)").
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         (i.e. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho_CT25   =  in-situ density                                [ kg m^-3 ]
%  in_funnel  =  0, if (SA, CT and p) are outside the "funnel" 
%             =  1, if (SA, CT and p) are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (23rd July, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See appendix A.20 and appendix K of this TEOS-10 Manual. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_rho_CT25:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_CT25: SA and CT must have same dimensions')
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
    error('gsw_rho_CT25: Inputs array dimensions arguments do not agree')
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

% These few lines ensure that SA is non-negative.
[I_neg_SA] = find(SA < 0);
if ~isempty(I_neg_SA)
    SA(I_neg_SA) = 0;
end

in_funnel = gsw_infunnel(SA,CT,p);

CT2 = CT.*CT; 

anum =             9.9984380290708214e+002 + ...
             CT.*( 7.1188090678940910e+000 + ...
             CT.*(-1.9459922513379687e-002 + ...
             CT *  6.1748404455874641e-004)) + ...
             SA.*( 2.8925731541277653e+000 + ...
             CT *  2.1471495493268324e-003 + ...
             SA *  1.9457531751183059e-003) + ...
              p.*( 1.1930681818531748e-002 + ...
            CT2 *  2.6969148011830758e-007 + ...
             SA *  5.9355685925035653e-006 + ...
              p.*(-2.5943389807429039e-008 + ...
            CT2 * -7.2734111712822707e-012));          
            
adenom =             1.00000000d+00 + ...
             CT.*( 7.0547681896071576e-003 + ...
             CT.*(-1.1753695605858647e-005 + ...
             CT.*( 5.9219809488274903e-007 + ...
             CT *  3.4887902228012519e-010))) + ...
             SA.*( 2.0777716085618458e-003 + ...
             CT.*(-2.2210857293722998e-008 + ...
            CT2 * -3.6628141067895282e-010) + ...
         sqrt(SA).*( 3.4688210757917340e-006 + ...
            CT2 *  8.0190541528070655e-010)) + ...
              p.*( 6.8314629554123324e-006 + ...
      (p.*CT).*(CT2 * -8.5294794834485446e-017 + ...
              p * -9.2275325145038070e-018));

rho_CT25 = anum./adenom;

%--------------------------------------------------------------------------
% This function calculates rho_CT25 using the computationally-efficient
% 25-term expression for density in terms of SA, CT and p. If one wanted to
% compute rho from SA, CT, and p with the full TEOS-10 Gibbs function,  
% the following lines of code will enable this.
%
%    pt0 = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t = gsw_pt_from_t(SA,pt0,pr0,p);
%    rho = gsw_rho(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%    rho = gsw_rho_CT(SA,CT,p)
%
%    or call the following, it is identical to the lines above.
%
%   [rho, ~, ~] = gsw_rho_alpha_beta_CT(SA,CT,p)
%
%-----------------This is the end of the alternative code------------------

if transposed
    rho_CT25 = rho_CT25';
end

end
