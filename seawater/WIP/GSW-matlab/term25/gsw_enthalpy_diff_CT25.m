function [enthalpy_diff_CT25, in_funnel] = gsw_enthalpy_diff_CT25(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff_CT25            difference of enthalpy at two pressures
%==========================================================================
%
% USAGE:
%  [enthalpy_diff, in_funnel] = gsw_enthalpy_diff_CT25(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between 
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  This 
%  function uses the computationally-efficient 25-term expression for
%  density in terms of SA, CT and p (McDougall et al., 2010).  The output
%  (enthalpy_diff) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow). 
%
% INPUT:
%  SA         =  Absolute Salinity                                 [ g/kg ]
%  CT         =  Conservative Temperature                         [ deg C ]
%  p_shallow  =  upper sea pressure                                [ dbar ]
%                (ie. shallower absolute pressure - 10.1325 dbar) 
%  p_deep     =  lower sea pressure                                [ dbar ]
%                (ie. deeper absolute pressure - 10.1325 dbar)
%
%  p_shallow and p_deep may have dimensions Mx1 or 1xN or MxN, 
%   where SA and CT are MxN.
%
% OUTPUT:
%  enthalpy_diff_CT25  =  difference of specific enthalpy          [ J/kg ]
%                         (deep minus shallow)
%  in_funnel           =  0, if SA, CT and p are outside the "funnel" 
%                      =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR: 
%   Trevor McDougall & Claire Roberts-Thomson.       [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) and (A.30.6) of this TEOS-10 Manual. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 4)
    error('gsw_enthalpy_diff_CT25: requires four inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (ms~=mt) | (ns~=nt)
    error('gsw_enthalpy_diff_CT25: SA & CT need to have the same dimensions')
end

if (mpu == 1) & (npu == 1)                           % p_shallowis a scalar 
    error('gsw_enthalpy_diff_CT25: need more than one pressure')
elseif (ns == npu) & (mpu == 1)                  % p_shallow is row vector,
    p_shallow = p_shallow(ones(1,ms), :);          % copy down each column.
elseif (ms == mpu) & (npu == 1)               % p_shallow is column vector,
    p_shallow = p_shallow(:,ones(1,ns));            % copy across each row.
elseif (ns == mpu) & (npu == 1)          % p_shallow is a transposed row vector,
    p_shallow = p_shallow';                              % transposed then
    p_shallow = p_shallow(ones(1,ms), :);                % copy down each column.
elseif (ms == mpu) & (ns == npu)
    % ok
end

if (mpl == 1) & (npl == 1)                             % p_deep is a scalar  
    error('gsw_enthalpy_diff_CT25: need more than one pressure')
elseif (ns == npl) & (mpl == 1)                     % p_deep is row vector,
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (npl == 1)                  % p_deep is column vector,
    p_deep = p_deep(:,ones(1,ns));                  % copy across each row.
elseif (ns == mpl) & (npl == 1)          % p_deep is a transposed row vector,
    p_deep = p_deep';                              % transposed then
    p_deep = p_deep(ones(1,ms), :);                % copy down each column.
elseif (ms == mpl) & (ns == npl)
    % ok
else
    error('gsw_enthalpy_diff_CT25: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA';
    CT = CT';
    p_shallow = p_shallow';
    p_deep = p_deep';
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

in_funnel_shallow = gsw_infunnel(SA,CT,p_shallow);
in_funnel_deep = gsw_infunnel(SA,CT,p_deep);
in_funnel = in_funnel_shallow.*in_funnel_deep;

db2Pa = 1e4; 

CT2 = CT.*CT; 
CT3 = CT.*CT2;

a0 = 1 + CT.*( 7.0547681896071576e-3 +...
         CT.*(-1.1753695605858647e-5 + ...
         CT.*(5.9219809488274903e-7 + ...
         CT.*3.4887902228012519e-10))) + ...
         SA.*( 2.0777716085618458e-3 + ...
         CT.*(-2.2210857293722998e-8 + ...
         CT2.*-3.6628141067895282e-10) + ...
    sqrt(SA).*(3.4688210757917340e-6 + ...
         CT2.*8.0190541528070655e-10));
a1 =  6.8314629554123324e-6;
a2 =  CT3*-8.5294794834485446e-17;
a3 =  CT*-9.2275325145038070e-18;

b0 =      9.9984380290708214e2 + ...
    CT.* (7.1188090678940910e0 + ...
    CT.*(-1.9459922513379687e-2 + ...
    CT.*  6.1748404455874641e-4)) + ...
    SA.*( 2.8925731541277653e0 + ...
    CT.*  2.1471495493268324e-3 + ...
    SA.*  1.9457531751183059e-3);
b1 = 0.5*(1.1930681818531748e-2 + ...
     CT2.*2.6969148011830758e-7 + ...
      SA.*5.9355685925035653e-6);
b2 = CT2.*-7.2734111712822707e-12 - 2.5943389807429039e-8;
b1sq = b1.*b1; 
sqrt_disc = sqrt(b1sq - b0.*b2);

N = a0 + (2*a3.*b0.*b1./b2 - a2.*b0)./b2;

M = a1 + (4*a3.*b1sq./b2 - (a3.*b0 + 2*a2.*b1))./b2;

A = b1 - sqrt_disc;
B = b1 + sqrt_disc;
delta_p = p_deep - p_shallow;
p_sum = p_deep + p_shallow;
part1 = b0 + p_shallow.*(2*b1 + b2.*p_shallow);

part2 = (B + b2.*p_deep).*(A + b2.*p_shallow);

part3 = (N.*b2 - M.*b1)./(b2.*(B - A));

enthalpy_diff_CT25 = db2Pa.*(delta_p.*(a2 - 2*a3.*b1./b2 + 0.5*a3.*p_sum)./b2 + ...
                     (M./(2*b2)).*log(1 + delta_p.*(2*b1 + b2.*p_sum)./part1) + ... 
                     part3.*log(1 + delta_p.*b2.*(B - A)./part2));

%--------------------------------------------------------------------------
% This function "gsw_enthalpy_diff_CT25.m" calculates enthalpy_diff_CT25 
% using the computationally-efficient 25-term expression for density in 
% terms of SA, CT and p.  If one wanted to compute the enthalpy difference 
% using the full TEOS-10 Gibbs function, the following lines of code will 
% enable this.
%
%    pt = gsw_pt_from_CT(SA,CT);
%    pr0 = zeros(size(SA)); 
%    t_shallow = gsw_pt_from_t(SA,pt,pr0,p_shallow);
%    t_deep = gsw_pt_from_t(SA,pt,pr0,p_deep);
%    enthalpy_diff = gsw_enthalpy(SA,t_deep,p_deep) - ...
%                    gsw_enthalpy(SA,t_shallow,p_shallow);
%
%    or call the following, it is identical to the lines above.
%
%    enthalpy_diff_CT = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep)
%
%-----------------This is the end of the alternative code------------------

  if transposed
      enthalpy_diff_CT25 = enthalpy_diff_CT25';
      in_funnel = in_funnel';
  end

end
