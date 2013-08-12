function [ outx, outy, doutx, douty ] = ConstR2( tx, ty )

outx = ones( size( tx ) );
outy = ones( size( ty ) );
doutx = zeros( size( tx ) );
douty = zeros( size( ty ) );