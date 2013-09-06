function [ outx, outy, doutx, douty ] = vdpR2( tx, ty )

scaling = 1.2;
vdpscale = 1;

sx = scaling * tx;
sy = scaling * ty;

outx = sy .* ( pi^2 - tx.^2 );
outy = -( sx + sy .* vdpscale .* ( sx.^2 - 1 ) ).* ( pi^2 - ty.^2 );

syms x y

dx = matlabFunction( diff( scaling * y * ( pi^2 - x^2 ), x ) );
dy = matlabFunction( diff( -( scaling * x + scaling * y * vdpscale ...
    * ( ( scaling * x )^2 - 1 ) ) * ( pi^2 - y^2 ), y ), 'vars', {x,y} );

doutx = dx( tx, ty );
douty = dy( tx, ty );