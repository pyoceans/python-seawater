function [SA_iref_cast, CT_iref_cast, p_iref_cast] = gsw_interp_McD_Klocker(spycnl,A)
% gsw_interp_McD_Klocker         linear interpolation of the reference cast
%==========================================================================
% This function interpolates the reference cast with respect to the 
% interpolating variable "spycnl".   This reference cast is at the location 
% 188E,4N from the reference data set which underlies the Jackett & 
% McDougall (1997) Neutral Density computer code.  This function finds the 
% values of SA, CT and p on this reference cast which correspond to the 
% value of isopycnal which is passed to this function from the function 
% "gsw_geo_str_McD_Klocker". The isopycnal could be either gamma_n or 
% sigma_2. if A is set to any of the following 's2','S2','sigma2','sigma_2'
% the interpolation will take place in sigma 2 space, any other input
% will result in the programme working in gamma_n space.
%
% REFERENCE:
%  Jackett, D. R. and T. J. McDougall, 1997: A neutral density variable
%   for the world’s oceans. Journal of Physical Oceanography, 27, 237-263.
%
% This fuction was adapted from Matlab's interp1q.
%==========================================================================

spycnl = spycnl(:);

if ~exist('A','var')
    A = 'gn';
elseif ~ischar(A)
    A = 'gn';
end %if
if strcmpi(A,'s2') == 1 | strcmpi(A,'s_2') == 1 | ...
        strcmpi(A,'sigma2') == 1 | strcmpi(A,'sigma_2') == 1 |...
        strcmpi(A,'s@') == 1 | strcmpi(A,'s_@') == 1 | ...
        strcmpi(A,'sigma@') == 1 | strcmpi(A,'sigma_@') == 1
    A = 's2';
end

temp = which('gsw_data_v2_0.mat');
dir = temp(1:((length(temp) - 17)));

if strcmpi(A,'s2')
    load ([dir,'gsw_data_v2_0'],'SA_ref_cast','CT_ref_cast','p_ref_cast','sigma_2_ref_cast');
    spycnl_ref_cast = sigma_2_ref_cast;
else
    load ([dir,'gsw_data_v2_0'],'SA_ref_cast','CT_ref_cast','p_ref_cast','gamma_n_ref_cast');
    spycnl_ref_cast = gamma_n_ref_cast;
end

[min_spycnl_ref_cast,Imin_spycnl_ref_cast] = min(spycnl_ref_cast);
[Ishallow] = find(spycnl <= min_spycnl_ref_cast); % Set equal to the shallowest bottle.
SA_iref_cast(Ishallow) = SA_ref_cast(Imin_spycnl_ref_cast);
CT_iref_cast(Ishallow) = CT_ref_cast(Imin_spycnl_ref_cast);
p_iref_cast(Ishallow) = p_ref_cast(Imin_spycnl_ref_cast);

[max_spycnl_ref_cast,Imax_spycnl_ref_cast] = max(spycnl_ref_cast);
[Ideep] = find(spycnl >= max_spycnl_ref_cast);       % Set equal to the deepest bottle.
SA_iref_cast(Ideep) = SA_ref_cast(Imax_spycnl_ref_cast);
CT_iref_cast(Ideep) = CT_ref_cast(Imax_spycnl_ref_cast);
p_iref_cast(Ideep) = p_ref_cast(Imax_spycnl_ref_cast);

[I] = find(spycnl >= 21.805 & spycnl <= 28.3614);

xi = spycnl(I);

x = spycnl_ref_cast;

siz = size(xi);
if ~isscalar(xi)
   [xxi, k] = sort(xi);
   [dum, j] = sort([x;xxi]);
   r(j) = 1:length(j);
   r = r(length(x)+1:end) - (1:length(xxi));
   r(k) = r;
   r(xi==x(end)) = length(x)-1;
   ind = find((r>0) & (r<length(x)));
   ind = ind(:);
   SA_ref_casti = NaN(length(xxi),size(SA_ref_cast,2),superiorfloat(x,SA_ref_cast,xi));
   CT_ref_casti = NaN(length(xxi),size(CT_ref_cast,2),superiorfloat(x,CT_ref_cast,xi));
   p_ref_casti = NaN(length(xxi),size(p_ref_cast,2),superiorfloat(x,p_ref_cast,xi));
   rind = r(ind);
   xrind = x(rind);
   u = (xi(ind)-xrind)./(x(rind+1)-xrind);
   SArind = SA_ref_cast(rind,:);
   CTrind = CT_ref_cast(rind,:);
   prind = p_ref_cast(rind,:);
   SA_ref_casti(ind,:) = SArind + bsxfun(@times,SA_ref_cast(rind+1,:)-SArind,u);
   CT_ref_casti(ind,:) = CTrind + bsxfun(@times,CT_ref_cast(rind+1,:)-CTrind,u);   
   p_ref_casti(ind,:) = prind + bsxfun(@times,p_ref_cast(rind+1,:)-prind,u);
else
   % Special scalar xi case
   r = find(x <= xi,1,'last');
   r(xi==x(end)) = length(x)-1;
   if isempty(r) || r<=0 || r>=length(x)
      SA_ref_casti = NaN(1,size(SA,2),superiorfloat(x,SA,xi));
      CT_ref_casti = NaN(1,size(CT,2),superiorfloat(x,CT,xi));
      p_ref_casti = NaN(1,size(p,2),superiorfloat(x,p,xi));
      
   else
      u = (xi-x(r))./(x(r+1)-x(r));
      SAr_ref_cast = SA_ref_cast(r,:);
      CTr_ref_cast = CT_ref_cast(r,:);
      pr_ref_cast = p_ref_cast(r,:);
      SA_ref_casti = SAr_ref_cast + bsxfun(@times,SA_ref_cast(r+1,:)-SAr_ref_cast,u);
      CT_ref_casti = CTr_ref_cast + bsxfun(@times,CT_ref_cast(r+1,:)-CTr_ref_cast,u);      
      p_ref_casti = pr_ref_cast + bsxfun(@times,p_ref_cast(r+1,:)-pr_ref_cast,u);
   end
end

if min(size(SA_ref_casti)) == 1 && numel(xi) > 1
   SA_ref_casti = reshape(SA_ref_casti,siz);
   CT_ref_casti = reshape(CT_ref_casti,siz);
   p_ref_casti = reshape(p_ref_casti,siz);
end

SA_iref_cast(I) = SA_ref_casti;
CT_iref_cast(I) = CT_ref_casti;
p_iref_cast(I) = p_ref_casti;

end