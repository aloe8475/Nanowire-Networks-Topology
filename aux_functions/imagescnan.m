function handle=imagescnan(IM)
% function imagescnan(IM)
% -- to display NaNs in imagesc as white/black

% white
nanjet = [ 1,1,1; jet  ];
nanjetLen = length(nanjet); 
pctDataSlotStart = 2/nanjetLen;
pctDataSlotEnd   = 1;
pctCmRange = pctDataSlotEnd - pctDataSlotStart;

dmin = nanmin(IM(:));
dmax = nanmax(IM(:));
dRange = dmax - dmin;   % data range, excluding NaN

cLimRange = dRange / pctCmRange;
cmin = dmin - (pctDataSlotStart * cLimRange);
cmax = dmax;
handle=imagesc(IM);
set(gcf,'colormap',nanjet);
caxis([cmin cmax]);
end 

function h = imlegend(h, dv, str, shapevalue)

% h is gcf
% dv is discrete value
% str is legend
% shapevalue can adjust the shape of the axes. which contain 4 values.
%   1st         movement in horz [no move is 0]
%   2nd         movement in vert [no move is 0]
%   3rd         zoom in horz [no zoom is 1]
%   4th         zoom in vert [no zoom is 1]
%
% imagesc legend
%

% window size
if nargin < 4
    shapevalue = [0, 0, 1, 1];
end

horzmove = shapevalue(1);
vertmove = shapevalue(2);
horzzoom = shapevalue(3);
vertzoom = shapevalue(4);

cdt = get(h, 'CData');
% otherix = ~multifind(cdt, dv);
% cdt(otherix) = nan;
cdtmin = min(cdt(:));
cdtmax = max(cdt(:));

% if interval of dv is larger than [cmin, cmax],
% than put min(dv) and max(dv) instead of cmin&cmax.
dvmin = min(dv);
dvmax = max(dv);

cmin = min(cdtmin, dvmin);
cmax = max(cdtmax, dvmax);

% axes position
ap1_l = 0.1 + horzmove;
ap1_b = 0.1 + vertmove;
ap1_len = 0.7 * horzzoom;
ap1_hei = 0.8 * vertzoom;

horzblank = 0.05 * horzzoom;
vertblank = 0.01 * vertzoom;
ap2_l = ap1_l + ap1_len + horzblank;
ap2_b = ap1_b;
ap2_len = 0.05 * horzzoom;
ap2_hei = 0.05 * vertzoom;

him = gca;
% [0.13 0.11 0.775 0.815]
% [0.8246 0.1143 0.0446 0.8167]
set(him, 'position', [ap1_l, ap1_b, ap1_len, ap1_hei])
set(him, 'Clim', [cmin, cmax])

n = length(dv);
hleg = nan(1, n);

for i = 1:n
    if dv(i) < cdtmin || dv(i) > cdtmax
        continue
    end
    hleg(i) = axes('position', [ap2_l, ap2_b + (i -1) * (ap2_hei + vertblank), ap2_len, ap2_hei]); %#ok<LAXES>
    imagesc(dv(i), [cmin, cmax])
    set(gca, 'xtick', [], 'ytick', [])
    xlim = get(gca, 'XLim');
    ylim = get(gca, 'YLim');
    text(xlim(2), (ylim(1) + ylim(2))/2, [' ', str{i}], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left')
end

axes(him) % set focus on axis.

end