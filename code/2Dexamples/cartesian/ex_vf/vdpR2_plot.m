function dy = vdpR2_plot( t, y )

scaling = 1.2;
vdpscale = 1;

dy = zeros( 2, 1 );

sx = scaling * y( 1 );
sy = scaling * y( 2 );

dy( 1 ) = sy * ( pi^2 - y( 1 ) );
dy( 2 ) = -( sx + sy * vdpscale * ( sx^2 - 1 ) )  * ( pi^2 - y( 2 ) );

