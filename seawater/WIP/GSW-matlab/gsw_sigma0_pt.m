function sigma0_pt = gsw_sigma0_pt(SA,pt0)

% gsw_sigma0                                     Potential density anomaly, 
%                                 being potential density minus 1000 kg/m^3 
% =========================================================================
%
% USAGE:
%   sigma0_pt = gsw_sigma0(SA,pt0)
%
% DESCRIPTION:
%  Calculates potential density anomaly with reference sea pressure of 
%  zero (0) dbar.  The temperature input to this function is potential
%  temperature referenced to zero dbar. 
%
% INPUT:
%  SA   =   Absolute Salinity                                      [ g/kg ]
%  pt0  =   potential temperature with respect to a               
%           reference sea pressure of 0 dbar (ITS-90)             [ deg C ]
%
%  SA & pt0 need to have the same dimensions.
%
% OUTPUT:
%  sigma0_pt  =  potential density anomaly with                     [ kg/m^3 ]
%                respect to a reference pressure of 0 dbar,   
%                that is, potential density minus 1000 kg/m^3.
%
% AUTHOR: 
%  Trevor McDougall & Paul Barker  [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0 (26th August, 2010)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.6.1) of this TEOS-10 Manual. 
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2) 
   error('gsw_sigma0: Requires two inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(pt0);

if (mt ~= ms | nt ~= ns)
    error('gsw_sigma0: SA and pt0 must have same dimensions')
end

if ms == 1
    SA = SA';
    pt0 = pt0';
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

sfac = 0.0248826675584615;                % sfac = 1/(40*(35.16504/35));
x2 = sfac.*SA; 
x  = sqrt(x2); 
y  = pt0.*0.025d0; 
   
g03 = 100015.695367145 + ...
       y.*(-270.983805184062 + ...
       y.*(1455.0364540468 + ...
       y.*(-672.50778314507 + ...
       y.*(397.968445406972 + ...
       y.*(-194.618310617595 + ...
       y.*(63.5113936641785 - ...
       y.*9.63108119393062))))));
                                                           
g08 = x2.*(-3310.49154044839 + ...
       x.*(199.459603073901 + ...
       x.*(-54.7919133532887 + ...
       x.* 36.0284195611086 - ...
       y.* 22.6683558512829) + ...
       y.*(-175.292041186547 + ...
       y.*(383.058066002476 + ...
       y.*(-460.319931801257 + ...
       y.* 234.565187611355)))) + ...
       y.*(729.116529735046 + ...
       y.*(-860.764303783977 + ...
       y.*(694.244814133268 + ...
       y.*(-297.728741987187)))));
     
sigma0_pt = 100000000./(g03 + g08) - 1000.0;

if transposed
   sigma0_pt = sigma0_pt';
end

% The output, being potential density anomaly, has units of kg/m^3 and is 
% this particular potential density (with referece pressure (pr) = 0 dbar) 
% and with 1000 kg/m^3 subtracted from it. 

end
