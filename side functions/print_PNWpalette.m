% Palette Print Function
function print_PNWpalette(x)
    pallength = size(x, 1);
    colormap(pnw_palette(x));
    colorbar('Location', 'southoutside');
    title(sprintf('%s, n=%d', x, pallength));
end