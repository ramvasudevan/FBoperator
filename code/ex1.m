NBases = 5;
compute_all_eigs = 1;
Xp = zeros( 2 * NBases + 1 );

for k = -NBases:NBases
    if ( k == -NBases )
        Xp( 2, k + NBases + 1 ) = 0.5 * ( k + 1 );
    elseif ( k == NBases )
        Xp( k + NBases, k + NBases + 1 ) = 0.5 * ( 1 - k );
    else
        Xp( k + NBases + 2, k + NBases + 1 ) = 0.5 * ( k + 1 );
        Xp( k + NBases, k + NBases + 1 ) = 0.5 * ( 1 - k );
    end
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
outp = zeros( 1, NSamples );
for j = 1:NSamples
    for k = -NBases:NBases
        %         if ( abs( V( k + NBases +1, zeroeig( 1 ) ) ) < 1e-10 )
        %             continue;
        %         end
        outp( j ) = outp( j ) + V( k + NBases + 1, zeroeig( 1 ) ) * exp( 1i * k * x( j ) );
    end
end

% figure; bar3( Xp );
figure; plot( x, outp );


