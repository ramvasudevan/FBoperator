function [ outx, outy, doutx, douty ] = SinR2( tx, ty )

outx = sin( tx );
outy = sin( ty );
doutx = cos( tx );
douty = cos( ty );