function [rho, alpha_wrt_CT, beta_const_CT, in_funnel] = gsw_rho_alpha_beta_CT25(SA,CT,p)

% gsw_rho_alpha_beta_CT25       In-situ density, thermal expansion & saline 
%                                contraction coefficient (25-term equation)
%==========================================================================
% 
% USAGE:  
% [rho, alpha_wrt_CT, beta_const_CT, in_funnel] = gsw_rho_alpha_beta_CT25(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from 
%  Absolute Salinity and Conservative Temperature.  This function uses the
%  computationally-efficient 25-term expression for density in terms of 
%  SA, CT and p (McDougall et al., 2010).
%
%  Note that potential density (pot_rho) with respect to reference pressure
%  pr is obtained by calling this function with the pressure argument being
%  pr as in [pot_rho, ~, ~, in_funnel] = gsw_rho_alpha_beta_CT25(SA,CT,pr).
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%          (ie. absolute pressure - 10.1325 dbar)
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  rho            =  in-situ density                            [ kg m^-3 ]
%  alpha_wrt_CT   =  thermal expansion coefficient                  [ 1/K ]
%                    with respect to Conservative Temperature
%  beta_const_CT  =  saline (i.e. haline) contraction coefficient  [ kg/g ]
%                    at constant Conservative Temperature
%  in_funnel      =  0, if SA, CT and p are outside the "funnel" 
%                 =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker   [ help_gsw@csiro.au ]
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
   error('gsw_rho_alpha_beta_CT25:  Requires three inputs')
end %if
if ~(nargout == 3 | nargout == 4)
   error('gsw_rho_alpha_beta_CT25:  Requires three or four outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_alpha_beta_CT25: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
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
    error('gsw_rho_alpha_beta_CT25: Inputs array dimensions arguments do not agree')
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
sqrtSA = sqrt(SA); 
pCT = p.*CT;

anum =             9.9984380290708214e+002 + ...
             CT.*( 7.1188090678940910e+000 + ...
             CT.*(-1.9459922513379687e-002 + ...
             CT *  6.1748404455874641e-004)) + ...
             SA.*( 2.8925731541277653e+000 + ...
             CT *  2.1471495493268324e-003 + ...
             SA *  1.9457531751183059e-003) + ...
              p.*( 1.1930681818531748e-002+ ...
            CT2 *  2.6969148011830758e-007 + ...
             SA *  5.9355685925035653e-006 + ...
              p.*(-2.5943389807429039e-008 + ...
            CT2 * -7.2734111712822707e-012));          
            
adenom =           1.00000000e+00 + ...
             CT.*( 7.0547681896071576e-003 + ...
             CT.*(-1.1753695605858647e-005 + ...
             CT.*( 5.9219809488274903e-007 + ...
             CT *  3.4887902228012519e-010))) + ...
             SA.*( 2.0777716085618458e-003 + ...
             CT.*(-2.2210857293722998e-008 + ...
            CT2 * -3.6628141067895282e-010) + ...
         sqrtSA.*( 3.4688210757917340e-006 + ...
            CT2.*  8.0190541528070655e-010))+ ...
              p.*( 6.8314629554123324e-006 + ...
      pCT.*(CT2.* -8.5294794834485446e-017 + ...
              p.* -9.2275325145038070e-018));

rec_adenom = ones(size(adenom))./adenom;

rho = anum.*rec_adenom;

rec_rho = ones(size(rho))./rho; 
   
anum_CT =          7.118809067894091e+00 + ...
             CT.*(-3.8919845026759374e-02 + ...
             CT.*  1.8524521336762394e-03) + ...
             SA.*  2.1471495493268324e-03 + ...
            pCT.*( 5.393829602366152e-07 - ...
              p.*  1.454682234256454e-11);
          
adenom_CT =        7.0547681896071576d-03 + ...
             CT.*(-2.35073912117172980d-05 + ...
             CT.*( 1.7765942846482467d-06 + ...
             CT.*  1.3955160891205007d-09)) + ...
             SA.*(-2.2210857293722998d-08 - ...
            CT2.*  1.09884423203685860d-09 + ...
     CT.*sqrtSA.*  1.6038108305614131d-09) - ...
     p.*p.*(CT2.*  2.5588438450345636d-16 + ...
              p.*  9.227532514503807d-18);
          
alpha_wrt_CT = (adenom_CT - anum_CT.*rec_rho).*rec_adenom;

anum_SA =          2.8925731541277653 + ...
             CT.*  2.1471495493268324e-03 + ...
             SA.*  3.8915063502366117e-03 + ...
              p.*  5.935568592503565e-06;

adenom_SA =        2.077771608561846e-03 + ...
             CT.*(-2.2210857293722998e-08 - ...
            CT2.*  3.6628141067895287e-10) + ...
         sqrtSA.*( 5.203231613687601e-06 + ...
            CT2.*  1.2028581229210597e-09);
          
beta_const_CT = (anum_SA.*rec_rho - adenom_SA).*rec_adenom;

%--------------------------------------------------------------------------
% This function calculates rho, alpha_wrt_CT and beta_const_CT using the 
% computationally-efficient 25-term expression for density in terms of SA, 
% CT and p. If one wanted to compute rho, alpha_wrt_CT and beta_const_CT 
% with the full TEOS-10 Gibbs function expression for density, the 
% following lines of code will do this.
%
%    pt0 = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t = gsw_pt_from_t(SA,pt0,pr0,p);
%    rho = gsw_rho(SA,t,p);
%    alpha_wrt_CT = gsw_alpha_wrt_CT(SA,t,p);
%    beta_const_CT = gsw_beta_const_CT(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%   [rho, alpha_wrt_CT, beta_const_CT] = gsw_rho_alpha_beta_CT(SA,CT,p)
%
%--------------This is the end of the alternative code---------------------

if transposed
    rho = rho';
    alpha_wrt_CT = alpha_wrt_CT';
    beta_const_CT = beta_const_CT';
    in_funnel = in_funnel';
end

end
