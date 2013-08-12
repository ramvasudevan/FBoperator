function out = sampvfS1( s )

bump = 0.5;

out = zeros( size( s ) );
idxpb = find( s > bump );
idxmb = find( s < -bump );
idxob = find( s <= bump & s >= -bump );

% out( idxpb ) = exp( -( s( idxpb ) - bump ).^4 ) .* sin( s( idxpb ) );
% out( idxmb ) = exp( -( s( idxmb ) + bump ).^4 ) .* sin( s( idxmb ) );
% out( idxob ) = sin( s( idxob ) );

out( idxpb ) = sin( s( idxpb ) ).^4 .* sin( s( idxpb ) + bump );
out( idxmb ) = sin( s( idxmb ) ).^4 .* sin( s( idxmb ) + bump );
out( idxob ) = sin( s( idxob ) ).^4 .* sin( s( idxob ) + bump );

if( find( s == -pi ) )
    out(  s == -pi  ) = 0;
elseif( find( s == pi ) )
    out( s == pi  ) = 0;
end