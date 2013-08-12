function out = deriv_sampvfS1( s )

bump = 0.5;

idxpb = find( s > bump );
idxmb = find( s < -bump );
idxob = find( s <= bump & s >= -bump );

syms x
dp = matlabFunction( diff( sin( x )^4 * sin( x + bump ), x ) );
dm = matlabFunction( diff( sin( x )^4 * sin( x + bump ), x ) );
do = matlabFunction( diff( sin( x )^4 * sin( x + bump ), x ) );

% out( idxpb ) = (4 .* exp( -( 1 - 2 .* atanh( th( idxpb )/pi ) ).^2 ) .* atanh( th( idxpb )/pi ).^2 .* ( 8 .* pi .* atanh( th( idxpb )/pi ).^2 + (2 .* th( idxpb ) - 4 .* pi ).* atanh( th( idxpb )/pi ) - 3 .* pi ) )/pi;
% out( idxmb ) = (4 .* exp( -( 2 .* atanh( th( idxmb )/pi ) + 1 ) .^2 ) .* atanh( th( idxmb )/pi ).^2 .* ( 8 .* pi .* atanh( th( idxmb )/pi ).^2 + 2 .* ( th( idxmb ) + 2 .* pi).* atanh( th( idxmb )/pi ) - 3 .* pi ) )/pi;
% out( idxob ) = (4 .* atanh( th( idxob )/pi ).^2 .* (2 .* th( idxob ) .* atanh( th( idxob )/pi ) - 3 .* pi ) )/pi;

% out( idxpb ) = -th( idxpb )/pi .* exp( -( s( idxpb ) - bump ).^2 ) .* -( s( idxpb ) + 1 ).^5 + exp( -( s( idxpb ) - bump ).^2 ) .* ( s( idxpb ) + 1 ).^4 .* ( 2 .* s( idxpb ).^2 - 8 .* s( idxpb ) - 15 );
% out( idxmb ) = -th( idxmb )/pi .* exp( -( s( idxmb ) + bump ).^2 ) .* -( s( idxmb ) + 1 ).^5 + exp( -( s( idxmb ) + bump ).^2 ) .* ( s( idxmb ) + 1 ).^4 .* ( 2 .* s( idxmb ).^2 + 12 .* s( idxmb ) + 5 );
% out( idxob ) = -th( idxob )/pi .* -( s( idxob ) + 1 ).^5 - 5.* ( s( idxob ) + 1 ).^4;

out( idxpb ) = dp( s( idxpb ) );
out( idxmb ) = dm( s( idxmb ) );
out( idxob ) = do( s( idxob ) );

if( find( s == -pi ) )
    out(  s == -pi  ) = 0;
elseif( find( s == pi ) )
    out(  s == pi  ) = 0;
end