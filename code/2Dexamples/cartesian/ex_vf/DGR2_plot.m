function dy = DGR2_plot( t, y )

dy = zeros( 2, 1 );

dy( 1 ) = -pi * sin( y( 1 ) ) * cos( y( 2 ) );
dy( 2 ) = pi * cos( y( 1 ) ) * sin( y( 2 ) );


