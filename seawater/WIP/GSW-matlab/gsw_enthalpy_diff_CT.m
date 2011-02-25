function enthalpy_diff_CT = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep)

% gsw_enthalpy_diff_CT              difference of enthalpy at two pressures
%==========================================================================
%
% USAGE:
%  enthalpy_diff_CT = gsw_enthalpy_diff_CT(SA,CT,p_shallow,p_deep)
%
% DESCRIPTION:
%  Calculates the difference of the specific enthalpy of seawater between
%  two different pressures, p_deep (the deeper pressure) and p_shallow
%  (the shallower pressure), at the same values of SA and CT.  The output
%  (enthalpy_diff) is the specific enthalpy evaluated at (SA,CT,p_deep)
%  minus the specific enthalpy at (SA,CT,p_shallow).
%
%  Note that this function uses the full Gibbs function.  There is an
%  alternative to calling this function, namely
%  gsw_enthalpy_diff_CT25(SA,CT,p), which uses the computationally
%  efficient 25-term expression for density in terms of SA, CT and p
%  (McDougall et al., (2010)).  For SA, CT and p values which fall inside
%  the oceanographic "funnel" (McDougall et al., 2010), this
%  computationally efficient (i. e. faster) 25-term version fits the
%  underlying laboratory density data almost as well as does the density
%  derived from the full TEOS-10 Gibbs function.
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
%  enthalpy_diff_CT  =  difference of specific enthalpy          [ J/kg ]
%                         (deep minus shallow)
%
% AUTHOR:
%   Trevor McDougall and Paul Barker       [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (29th September, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqns. (3.32.2) of this TEOS-10 Manual.
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
    error('gsw_enthalpy_diff_CT: requires four inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mpu,npu] = size(p_shallow);
[mpl,npl] = size(p_deep);

if (ms~=mt) | (ns~=nt)
    error('gsw_enthalpy_diff_CT: SA & CT need to have the same dimensions')
end

if (mpu == 1) & (npu == 1)                           % p_shallowis a scalar
    error('gsw_enthalpy_diff_CT: need more than one pressure')
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
    error('gsw_enthalpy_diff_CT: need more than one pressure')
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
    error('gsw_enthalpy_diff_CT: Inputs array dimensions arguments do not agree')
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

pt = gsw_pt_from_CT(SA,CT);
pr0 = zeros(size(SA));
t_shallow = gsw_pt_from_t(SA,pt,pr0,p_shallow);
t_deep = gsw_pt_from_t(SA,pt,pr0,p_deep);
enthalpy_diff_CT = gsw_enthalpy(SA,t_deep,p_deep) - ...
                     gsw_enthalpy(SA,t_shallow,p_shallow);

if transposed
    enthalpy_diff_CT = enthalpy_diff_CT';
end

end
