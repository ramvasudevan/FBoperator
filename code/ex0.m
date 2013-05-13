NBases = 22;
compute_all_eigs = 1;
Xp = zeros( 2 * NBases + 1 );

for k = -NBases:NBases
    Xp( k + NBases + 1, k + NBases + 1 ) = -1i * k;
end

if ( compute_all_eigs )
    % the following code does not do a great job of computing the eigenvalues
    % probably due to some scaling issue, but the true reason is unclear
    % [ V, D ] = eig( Xp );
    % this works a bit better
    [ V, D ] = eig( Xp, eye( size( Xp ) ));
    
    zeroeig = [];
    for k = 1:length(D)
        if ( D( k, k ) == 0 )
            zeroeig = [ zeroeig; k ];
        end
    end
else
    % this works fast
    [ V, D ] = eigs( Xp, eye( size( Xp ) ), 1, 1e-5 );
    
    zeroeig = 1;
end

% zeroeig = 2 * NBases + 1;
NSamples = 1001;
x = linspace( -pi, pi, NSamples );
out = zeros( 1, NSamples );
for j = 1:NSamples
    for k = -NBases:NBases
        out( j ) = out( j ) + V( k + NBases + 1, zeroeig( 1 ) ) * exp( 1i * k * x( j ) );
    end
end

figure; plot( x, out );


