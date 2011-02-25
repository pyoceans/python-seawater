function [geo_strf_dyn_height, in_funnel] = gsw_geo_strf_dyn_height(SA,CT,p,delta_p,interp_style)

% gsw_geo_strf_dyn_height                            dynamic height anomaly
%==========================================================================
%
% USAGE:  
%  [geo_strf_dyn_height, in_funnel] = gsw_geo_strf_dyn_height(SA,CT,p,delta_p,interp_style)
%
% DESCRIPTION:
%  Calculates dynamic height anomaly as the integral of specific volume 
%  anomaly from the the sea surface pressure (0 Pa) to the pressure p.
%  Hence, geo_strf_dyn_height is the dynamic height anomaly with respect to
%  the sea surface.  This is the geostrophic streamfunction for the 
%  difference between the horizontal velocity at the pressure concerned, p,
%  and the horizontal velocity at the sea surface.  Dynamic height anomaly
%  is the geostrophic streamfunction in an isobaric surface.  The reference
%  values used for the specific volume anomaly are SA = SSO = 35.16504 g/kg
%  and CT = 0 deg C.  This function calculates specific volume anomaly 
%  using the computationally efficient 25-term expression for specific 
%  volume of McDougall et al. (2010). 
%  Under the default setting, this function evaluates the pressure integral
%  of specific volume using SA and CT “interploted” with respect to pressure
%  using a scheme based on the method of Reiniger and Ross (1968).  Our 
%  method uses a weighted mean of (i) values obtained from linear 
%  interpolation of the two nearest data points, and (ii) a linear 
%  extrapolation of the pairs of data above and below. This "curve fitting"
%  method resembles the use of cubic splines.  If the option “linear” is 
%  chosen, the function interpolates Absolute Salinity and Conservative 
%  Temperature linearly with presure in the vertical between “bottles”.
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  CT   =  Conservative Temperature                               [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%         ( ie. absolute pressure - 10.1325 dbar )
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OPTIONAL:
%  delta_p      = maximum interpolation distance between bottles.
%  interp_style = interpolation technique.
%               = if nothing is entered the programme defaults to "curved" 
%                 interpolation between bottles in the vertical.
%               = if "linear" or "lin" is entered then the programme 
%                 interpolates linearly between bottles in the
%                 vertical.
%                
% OUTPUT:
%  geo_strf_dyn_height = dynamic height anomaly                 [ m^2/s^2 ]
%
%  in_funnel          =  0, if SA, CT and p are outside the "funnel" 
%                     =  1, if SA, CT and p are inside the "funnel"
%  Note. The term "funnel" describes the range of SA, CT and p over which 
%    the error in the fit of the computationally-efficient 25-term 
%    expression for density was calculated (McDougall et al., 2010).
%
% AUTHOR:  
%  Paul Barker, Jeff Dunn and Trevor McDougall [ help_gsw@csiro.au ]
%
% VERSION NUMBER: 2.0, 25th August, 2010 
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) and section 3.27 of this TEOS-10 Manual. 
%
%  McDougall T. J., D. R. Jackett, P. M. Barker, C. Roberts-Thomson, R.
%   Feistel and R. W. Hallberg, 2010:  A computationally efficient 25-term 
%   expression for the density of seawater in terms of Conservative 
%   Temperature, and related properties of seawater.  To be submitted 
%   to Ocean Science Discussions. 
%
%  Reiniger, R. F. and C. K. Ross, 1968: A method of interpolation with
%   application to oceanographic data. Deep-Sea Res. 15, 185-193.
% 
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 | nargin == 4 | nargin == 5)
   error('gsw_geo_strf_dyn_height:  Requires three, four or five inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (ms~=mt) | (ns~=nt)
    error('gsw_geo_strf_dyn_height: SA & CT need to have the same dimensions')
elseif (ms*ns == 1)
    error('gsw_geo_strf_dyn_height: There must be at least 2 values')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    error('gsw_geo_strf_dyn_height: need more than one pressure')
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
    error('gsw_geo_strf_dyn_height: Inputs array dimensions arguments do not agree')
end %if

if ~exist('delta_p','var')
    delta_p = 1;
end

if ~exist('interp_style','var')
    interp_style = 'curve';
elseif strcmpi('interp_style','linear') == 1 | strcmpi('interp_style','lin') == 1 |...
        strcmpi('interp_style','linaer') == 1 | strcmpi('interp_style','lnear') == 1
    interp_style = 'linear';
end

if ms == 1
    SA = SA';
    CT = CT';
    p = p';
    transposed = 1;
else
    transposed = 0;
end
[mp,np] = size(p);

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

in_funnel = gsw_infunnel(SA,CT,p);

db2Pa = 1e4;

Ishallow = 1:(mp-1);
Ideep = 2:mp;
d_p = (p(Ideep,:)-p(Ishallow,:));

if ~isempty(find(d_p <= 0))
    error('gsw_geo_strf_dyn_height: pressure must be monotonic')
end

geo_strf_dyn_height = nan(size(SA));

[Ibg] = find(d_p > delta_p);
[Inz] = find(p(1,:) ~=0);

if isempty(Ibg) & isempty(Inz)
    % resolution is ok & profile extends to the surface.
    B = gsw_specvol_CT25(SA,CT,p);
    C = gsw_enthalpy_SSO_0_CT25(p);
    
%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly using the computationally
% efficient 25-term expression for density in terms of SA, CT and p. If one
% wanted to compute dynamic height anomaly with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable 
% this.
%
%    B = gsw_specvol_CT(SA,CT,p);
%    SA_SO = 35.16504*ones(size(SA));
%    CT_0 = zeros(size(CT));
%    C = gsw_enthalpy_CT(SA_SO,CT_0,p);
%
% Further down the page is a second section which also needs to be
% activated in order to compute dynamic height anomaly with the full 
% TEOS-10 Gibbs function expression for density.
%---------------This is the end of the alternative code--------------------

    B_av = zeros(size(SA));
    B_av(2:mp,:) = 0.5*(B(1:(end-1),:) + B(2:end,:));    
    dp = zeros(size(SA));
    dp(2:mp,:) = d_p;
    D = B_av.*dp.*db2Pa;
    geo_strf_dyn_height = C - cumsum(D);    
else
    % will need to interplolate profiles, but will check each profile.
    for Iprofile = 1:np
        [Inn] = find(~isnan(p(:,Iprofile)));
        [Ibg_i] = find(d_p(:,Iprofile) > delta_p);
        if isempty(Ibg_i)
            if min(p(Inn,Iprofile)) > 0
                SA_i = SA(Inn(1),Iprofile);
                SA_i(2:length(Inn)+1) = SA(Inn,Iprofile);                
                CT_i = CT(Inn(1),Iprofile);
                CT_i(2:length(Inn)+1) = CT(Inn,Iprofile);
                p_i = 0;
                p_i(2:length(Inn)+1) = p(Inn,Iprofile);
                [dummy Iidata Ibdata] = intersect(p_i(2:end),p(:,Iprofile));
            else
                SA_i = SA(Inn,Iprofile);
                CT_i = CT(Inn,Iprofile);
                p_i = p(Inn,Iprofile);
                [dummy Iidata Ibdata] = intersect(p_i,p(:,Iprofile));
            end
        else
            p_i = nan(2*round(max(p(Inn,Iprofile)/delta_p)),1);
            if p(Inn,Iprofile) > 0
                p_i(1) = 0;
                p_i(2) = min(p(Inn,Iprofile));
                top_pad = 1;
                p_cnt = 2;
            else
                p_i(1) = min(p(Inn,Iprofile));
                top_pad = 0;
                p_cnt = 1;
            end
            for Ibottle = 1:(length(Inn)-1)
                dp_iIbottle = p(Inn(Ibottle+1),Iprofile) - p(Inn(Ibottle),Iprofile);
                p_iIbottle = p(Inn(Ibottle),Iprofile):dp_iIbottle/ceil(dp_iIbottle/delta_p):p(Inn(Ibottle+1),Iprofile);
                p_cnt_ld = p_cnt+1;
                p_cnt = p_cnt + length(p_iIbottle(2:length(p_iIbottle)));
                p_i(p_cnt_ld:p_cnt) = p_iIbottle(2:length(p_iIbottle));
            end
            p_i(p_cnt+1:end) = []; 
            p_i = p_i(:);

            if top_pad ==1
                Intrp = 2:length(p_i);
                SA_i = SA(Inn,Iprofile);
                CT_i = CT(Inn,Iprofile);
            else
                Intrp = 1:length(p_i);
            end
            [dummy Iidata Ibdata] = intersect(p_i(Intrp),p(:,Iprofile));

            if ~exist('interp_style','var')
                SA_i(Intrp) = pinterp_from_p(p(:,Iprofile),SA(:,Iprofile),p_i(Intrp));
                CT_i(Intrp) = pinterp_from_p(p(:,Iprofile),CT(:,Iprofile),p_i(Intrp));
                [Inan] = find(isnan(SA_i));
                if ~isempty(Inan)
                    [SA_i(Inan), CT_i(Inan)] = gsw_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan));
                end       
            elseif exist('interp_style','var') & strcmpi(interp_style,'linear')
                [SA_i(Intrp), CT_i(Intrp)] = gsw_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Intrp));
            else
                SA_i(Intrp) = pinterp_from_p(p(:,Iprofile),SA(:,Iprofile),p_i(Intrp));
                CT_i(Intrp) = pinterp_from_p(p(:,Iprofile),CT(:,Iprofile),p_i(Intrp));
                [Inan] = find(isnan(SA_i));
                if ~isempty(Inan)
                    [SA_i(Inan), CT_i(Inan)] = gsw_interp_SA_CT(SA(:,Iprofile),CT(:,Iprofile),p(:,Iprofile),p_i(Inan));
                end       
             end
        end
        B_i = gsw_specvol_CT25(SA_i(:),CT_i(:),p_i(:));
        C_i = gsw_enthalpy_SSO_0_CT25(p_i(2:end));

%--------------------------------------------------------------------------
% This function calculates dynamic height anomaly using the computationally
% efficient 25-term expression for density in terms of SA, CT and p. If one
% wanted to compute dynamic height anomaly with the full TEOS-10 Gibbs 
% function expression for density, the following lines of code will enable 
% this.
%
%    B_i = gsw_specvol_CT(SA_i,CT_i,p_i);
%    B_i = B_i(:);
%    SA_SO = 35.16504*ones(size(SA_i));
%    CT_0 = zeros(size(CT_i));
%    C_i = gsw_enthalpy_CT(SA_SO(2:end),CT_0(2:end),p_i(2:end));
%    C_i = C_i(:);
%
%---------------This is the end of the alternative code--------------------

        B_i_av = 0.5*(B_i(1:(end-1)) + B_i(2:end));
        Da_i = (B_i_av.*diff(p_i).*db2Pa);
        D_i(2:length(C_i)+1) = C_i - cumsum(Da_i);
        geo_strf_dyn_height(Ibdata,Iprofile) = D_i(Intrp(Iidata));
        %keyboard
        clear SA_i CT_i p_i 
    end
end

if transposed
   geo_strf_dyn_height = geo_strf_dyn_height';
   in_funnel = in_funnel';
end %if

end

%##########################################################################

function [sdat] = pinterp_from_p(odep,obs,sdep)
% pinterp_from_p.
%==========================================================================
% Interpolate values on arbitrary pressures (Designed for bottle data, but 
% fine for 2db CTD data because it handles any gaps safely).
% INPUT:
%  odep  - vector of pressures of the data.
%  obs   - corresponding data values, with nan indicating any gaps.
%  sdep  - pressure to interpolate to.
% OUTPUT:
%  sdat - interpolated values on at the required pressures.
% Jeff Dunn.   Copyright CSIRO Division of Marine Research.
%==========================================================================

global rr_int_cnt lin_int_cnt dir_sub_cnt r_extrp_cnt;
grad_lim = [];
maxdis = rr_int([],[],sdep);
odep = odep(:);
obs = obs(:);
sdep = sdep(:);
nlvl = length(sdep);
xfn = [0 300 1200 8000];
yfn = [7 15 75 150];
near_lim = interp1(xfn,yfn,sdep);
far_lim = 2*near_lim;
dir_lim = near_lim/5;
sdat = repmat(NaN,nlvl,1);
jj = find(isnan(obs) | isnan(odep) | odep<0 | odep>8000);
if ~isempty(jj)
    obs(jj) = [];
    odep(jj) = [];
end
if ~isempty(obs)
    jj = find((odep(2:end)-odep(1:end-1))<=0);
    if ~isempty(jj)
        obs(jj+1) = [];
        odep(jj+1) = [];
    end
end
ndeps = length(obs);
if ndeps == 0
    % RETURN if no data
    return;
end

if nargin<3 | isempty(maxdis); maxdis = 1; end

if ndeps < 4 | maxdis == -1
    sidx = (1:nlvl)';
else
    % RR INTERPOLATION
    sdat = rr_int(obs,odep,sdep,1,maxdis);
    sidx = find(isnan(sdat));
    rr_int_cnt = rr_int_cnt + nlvl - length(sidx);
end

if ~isempty(sidx)  & ndeps >= 2
    idx = sidx(find(sdep(sidx)>odep(1) & sdep(sidx)<odep(ndeps)));
    if ~isempty(idx)
        oidx = interp1(odep,1:ndeps,sdep(idx));
        dists = [sdep(idx)-odep(floor(oidx)) odep(ceil(oidx))-sdep(idx)];
        near = min(dists')';
        far = max(dists')';
        interp = idx(find(near<near_lim(idx) | far<far_lim(idx)));
        if ~isempty(interp)
            sdat(interp) = interp1(odep,obs,sdep(interp));
            sidx = find(isnan(sdat));
            lin_int_cnt = lin_int_cnt + length(interp);
        end
    end
end

if ~isempty(sidx)
    idx = round(interp1([-99999; odep; 99999],0:ndeps+1,sdep(sidx)));
    kk = find(abs(odep(idx)-sdep(sidx)) < near_lim(sidx));
    for jj = kk(:)'
        sdj = sdep(sidx(jj));
        odj = odep(idx(jj));
        x = sdj-odj;
        new = nan;
        if abs(x) > 1.5
            jll = find(abs(odep-sdj) < far_lim(sidx(jj)));
            if x > 0
                jll = flipud(jll);
            end
            if length(jll)<2 | max(abs(odep(jll)-odj)) < abs(x)
                jll = [];
            elseif any(abs(diff(odep(jll))) < 1.5)
                ll = jll(1);
                for mm = jll(2:end)'
                    if abs(odep(ll(end))-odep(mm))>1.5  & ...
                            (length(ll) < 4 | abs(odj-odep(mm)) < abs(x))
                        ll = [ll mm];
                    end
                end
                jll = ll;
            end
            if length(jll) >= 2
                if abs(max(obs(jll))-min(obs(jll)))<.005
                    new = obs(jll(1));
                else
                    xog = min(odep(jll));
                    cc = ([ones([length(jll) 1]) odep(jll)-xog]) \ obs(jll);
                    new = cc(1) + (sdj-xog)*cc(2);
                end
                r_extrp_cnt = r_extrp_cnt + 1;
                if ~isempty(grad_lim)
                    ofset = abs(obs(idx(jj))-new);
                    if ofset>abs(x*grad_lim) | ofset>offlim
                        new = nan;
                    end
                end
                sdat(sidx(jj)) = new;
            end
        end
        if isnan(new) & abs(x)<dir_lim(sidx(jj))
            sdat(sidx(jj)) = obs(idx(jj));
            dir_sub_cnt = dir_sub_cnt + 1;
        end
    end
    sidx = find(isnan(sdat));
end

end

%##########################################################################

function [yi,maxdis] = rr_int(y,x,xi,limchk,maxdis)
%==========================================================================
% References:
%  Reiniger RF & Ross CK. 1968.  A method for interpolation with application
%   to oceanographic data.  Deep-Sea Res. 15: 185-193
% Jeff Dunn  12/5/97  Copyright CSIRO Division of Marine Research
%==========================================================================

cfrac = 1/15;
coincid_frac = 1/200;
if nargin<4
    limchk = 1;
elseif isempty(limchk)
    limchk = 1;
end
limchk = 0;
xfn = [0 300 1800 8000];
yfn = [50 200 650 1250];
maxdis = interp1(xfn,yfn,xi);
maxdis = [maxdis(:) maxdis(:).*3];
nobs = length(x);
if nobs < 4
    yi = [];
    return;
end
if size(x,1) == 1;  x =x'; end
if size(y,1) == 1;  y =y'; end
if size(xi,1) == 1;  xi =xi'; end
tidx = (1:length(xi))';
yi = repmat(NaN,size(tidx));
if x(nobs)>x(1)
    tidx = tidx(find(xi>=x(2) & xi<=x(nobs-1)));
else
    tidx = tidx(find(xi>=x(nobs-1) & xi<=x(2)));
end
if ~isempty(tidx)
    oidx = interp1(x,(1:length(x))',xi(tidx));
    if ~isempty(tidx)
        ridx = round(oidx);
        coincid = find(abs(xi(tidx)-x(ridx)) < (min(maxdis(:))*coincid_frac));
        if ~isempty(coincid)
            yi(tidx(coincid)) = y(ridx(coincid));
            oidx(coincid) = [];
            tidx(coincid) = [];
        end
    end
    rej = find(oidx<2 | oidx>(length(x)-1));
    if ~isempty(rej)
        oidx(rej) = [];
        tidx(rej) = [];
    end
    if ~isempty(tidx)
        fidx = floor(oidx);
        cidx = ceil(oidx);
        inn = [xi(tidx)-x(fidx) x(cidx)-xi(tidx)];
        out = [x(fidx)-x(fidx-1) x(cidx+1)-x(cidx)];
        sumin = inn(:,1) + inn(:,2);
        outtest = abs(out(:,1) + out(:,2) + sumin);
        minsep = abs(min(out')'./sumin);
        rej = find(sumin>maxdis(tidx,1) | outtest>maxdis(tidx,2) | minsep<cfrac);
        if ~isempty(rej)
            tidx(rej) = [];
            inn(rej,:) = [];
            out(rej,:) = [];
            sumin(rej) = [];
            fidx(rej) = [];
            cidx(rej) = [];
        end
    end
    if ~isempty(tidx)
        % Calculate the inner interpolated and 2 outer extrapolated values:
        intp = y(fidx) + ((y(cidx)-y(fidx)).*inn(:,1)./sumin);
        ext1 = y(fidx) + ((y(fidx)-y(fidx-1)).*inn(:,1)./out(:,1));
        ext2 = y(cidx) + ((y(cidx)-y(cidx+1)).*inn(:,2)./out(:,2));
        % Construct the Reiniger&Ross reference curve equation
        % m = the power variable
        m = 1.7;
        top = (abs(ext1-intp).^m).*ext2 + (abs(intp-ext2).^m).*ext1;
        bot = abs(ext1-intp).^m + abs(intp-ext2).^m;
        kk = find(abs(bot)<1E-4);
        if ~isempty(kk)
            yi(tidx(kk)) = intp(kk);
            kk = find(abs(bot)>=1E-4);
            yi(tidx(kk)) = (intp(kk)+(top(kk)./bot(kk)))./2;
        else
            yi(tidx) = (intp+(top./bot))./2;
        end
        gt = y(fidx) > y(cidx);
        yi(tidx) = max([yi(tidx)'; y(fidx+gt)']);
        yi(tidx) = min([yi(tidx)'; y(cidx-gt)']);
    end
end
end

%##########################################################################

