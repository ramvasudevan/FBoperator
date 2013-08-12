function dy = vdpS1_plot( t, y )

scaling = 1.5;
vdpscale = 1;

dy = zeros( 2, 1 );

sx = 2 * atanh( scaling * y( 1 )/pi );
sy = 2 * atanh( scaling * y( 2 )/pi );

dy( 1 ) = ( pi^2 - y( 1 )^2 )/( 2 * pi ) * sy * ( pi^2 - y( 1 )^2 )^3 * ( pi^2 - y( 2 )^2 )^3;
dy( 2 ) = ( pi^2 - y( 2 )^2 )/( 2 * pi ) *-( sx + sy * vdpscale * ( sx^2 - 1) ) * ( pi^2 - y( 1 )^2 )^3 * ( pi^2 - y( 2 ).^2 )^3;



if( y( 1 ) == pi || y( 2 ) == pi )
    dy( 1 ) = 0;
    dy( 2 ) = 0;
elseif( y( 1 ) == -pi || y( 2 ) == -pi )
    dy( 1 ) = 0;
    dy( 2 ) = 0;
end
