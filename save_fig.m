function save_fig(h,fname)
%function save_fig(h,fname)
% wrapper for export_fig
%
% assumes plot is still up
% makes plot with a transparent background


set(h,'Color','w');
export_fig(h,fname); % matlab central function


end