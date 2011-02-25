function [rho, alpha_wrt_CT, beta_const_CT] = gsw_rho_alpha_beta_CT(SA,CT,p)

% gsw_rho_alpha_beta_CT25 (?)               in-situ density, thermal expansion
%                                  & saline contraction coefficient from CT
%==========================================================================
%
% USAGE:
% [rho, alpha_wrt_CT, beta_const_CT] = gsw_rho_alpha_beta_CT(SA,CT,p)
%
% DESCRIPTION:
%  Calculates in-situ density, the appropiate thermal expansion coefficient
%  and the appropriate saline contraction coefficient of seawater from
%  Absolute Salinity and Conservative Temperature.
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_rho_alpha_beta_CT25(SA,CT,p), which uses the computationally
%  efficient 25-term expression for density in terms of SA, CT and p
%  (McDougall et al., (2010)).  For SA, CT and p values which fall inside
%  the oceanographic "funnel" (McDougall et al., 2010), this
%  computationally efficient (i. e. faster) 25-term version fits the
%  underlying laboratory density data almost as well as does the density
%  derived from the full TEOS-10 Gibbs function.
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
%  beta_const_CT  =  saline contraction coefficient                [ kg/g ]
%                    at constant Conservative Temperature
%
% AUTHOR:
%  Trevor McDougall and Paul Barker           [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (30th September, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%     See sections (2.8), (2.18) and (2.19) of this TEOS-10 Manual.
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
   error('gsw_rho_alpha_beta_CT:  Requires three inputs')
end %if
if ~(nargout == 3)
   error('gsw_rho_alpha_beta_CT:  Requires three outputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_rho_alpha_beta_CT: SA and CT must have same dimensions')
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
    error('gsw_rho_alpha_beta_CT: Inputs array dimensions arguments do not agree')
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

pt0 = gsw_pt_from_CT(SA,CT);
pr0 = zeros(size(SA));
t = gsw_pt_from_t(SA,pt0,pr0,p);
rho = gsw_rho(SA,t,p);
alpha_wrt_CT = gsw_alpha_wrt_CT(SA,t,p);
beta_const_CT = gsw_beta_const_CT(SA,t,p);

if transposed
    rho = rho';
    alpha_wrt_CT = alpha_wrt_CT';
    beta_const_CT = beta_const_CT';
end

end
