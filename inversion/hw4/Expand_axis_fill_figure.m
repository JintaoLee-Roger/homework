function [ ] = Expand_axis_fill_figure( axis_handle )
% input:  the axis_handle of the target figure
% output: None
%
% Example: 
%       Expand_axis_fill_figure(gca)

% TightInset position
inset_vectior = get(axis_handle, 'TightInset');
inset_x = inset_vectior(1);
inset_y = inset_vectior(2);
inset_w = inset_vectior(3);
inset_h = inset_vectior(4);

% OuterPosition position
outer_vector = get(axis_handle, 'OuterPosition');
% Move Position origin to TightInset origin
pos_new_x = outer_vector(1) + inset_x; 
pos_new_y = outer_vector(2) + inset_y;
pos_new_w = outer_vector(3) - inset_w - inset_x; % reset Position's wide
pos_new_h = outer_vector(4) - inset_h - inset_y; % reset Position's height

% reset Position
set(axis_handle, 'Position', [pos_new_x, pos_new_y, pos_new_w, pos_new_h]);