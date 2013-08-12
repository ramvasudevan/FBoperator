function [ outx, outy, doutx, douty ] = DGR2( tx, ty )

outx = -pi * sin( tx ) .* cos( ty );
outy = pi * cos( tx ) .* sin( ty );
doutx = -pi * cos( tx ) .* cos( ty );
douty = pi * cos( tx ) .* cos( ty );