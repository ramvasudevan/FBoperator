function [ outx, outy, doutx, douty ] = vdpS1( tx, ty )

scaling = 1.5;
vdpscale = 1;

sx = 2 * atanh( scaling * tx/pi );
sy = 2 * atanh( scaling * ty/pi );

outx = ( pi^2 - tx.^2 )/( 2 * pi ) .* sy .* ( pi^2 - tx.^2 ).^2 .* ( pi^2 - ty.^2 ).^2;
outy = ( pi^2 - ty.^2 )/( 2 * pi ) .*-( sx + sy .* vdpscale .* ( sx.^2 - 1) ) .* ( pi^2 - tx.^2 ).^2 .* ( pi^2 - ty.^2 ).^2;

syms x y

dx = matlabFunction( diff( ( pi^2 - x^2 )/( 2 * pi ) * 2 * atanh( scaling * y/pi ) * ( pi^2 - x^2 )^2 * ( pi^2 - y^2 )^2, x ) );
dy = matlabFunction( diff( ( pi^2 - y^2 )/( 2 * pi ) *-( 2 * atanh( scaling * x/pi ) + 2 * atanh( scaling * y/pi ) * vdpscale ...
    * ( (2 * atanh( scaling * x/pi ) )^2 - 1 ) ) * ( pi^2 - x^2 )^2 * ( pi^2 - y^2 )^2, y ) );

doutx = dx( tx, ty );
douty = dy( tx, ty );

idxpi = tx == pi | ty == pi | tx == -pi | ty == -pi;

outx( idxpi ) = 0;
outy( idxpi ) = 0;
doutx( idxpi ) = 0;
douty( idxpi ) = 0;