function out = deriv_sampvf( th )

bump = 5;

out = zeros( size( th ) );
% s = 2 * atanh( th/pi );
% idxpb = find( s > bump );
% idxmb = find( s < -bump );
% idxob = find( s <= bump & s >= -bump );
% 
% syms x
% dp = matlabFunction( diff( ( pi^2 - x^2 )/( 2 * pi ) * exp( -( 2 * atanh( x/pi ) - bump )^2 ) * -( 2 * atanh( x/pi ) + 0.5 )^5, x ) );
% dm = matlabFunction( diff( ( pi^2 - x^2 )/( 2 * pi ) * exp( -( 2 * atanh( x/pi ) + bump )^2 ) * -( 2 * atanh( x/pi ) + 0.5 )^5, x ) );
% do = matlabFunction( diff( ( pi^2 - x^2 )/( 2 * pi ) * -( 2 * atanh( x/pi ) + 0.5 )^5, x ) );
% 
% % out( idxpb ) = (4 .* exp( -( 1 - 2 .* atanh( th( idxpb )/pi ) ).^2 ) .* atanh( th( idxpb )/pi ).^2 .* ( 8 .* pi .* atanh( th( idxpb )/pi ).^2 + (2 .* th( idxpb ) - 4 .* pi ).* atanh( th( idxpb )/pi ) - 3 .* pi ) )/pi;
% % out( idxmb ) = (4 .* exp( -( 2 .* atanh( th( idxmb )/pi ) + 1 ) .^2 ) .* atanh( th( idxmb )/pi ).^2 .* ( 8 .* pi .* atanh( th( idxmb )/pi ).^2 + 2 .* ( th( idxmb ) + 2 .* pi).* atanh( th( idxmb )/pi ) - 3 .* pi ) )/pi;
% % out( idxob ) = (4 .* atanh( th( idxob )/pi ).^2 .* (2 .* th( idxob ) .* atanh( th( idxob )/pi ) - 3 .* pi ) )/pi;
% 
% % out( idxpb ) = -th( idxpb )/pi .* exp( -( s( idxpb ) - bump ).^2 ) .* -( s( idxpb ) + 1 ).^5 + exp( -( s( idxpb ) - bump ).^2 ) .* ( s( idxpb ) + 1 ).^4 .* ( 2 .* s( idxpb ).^2 - 8 .* s( idxpb ) - 15 );
% % out( idxmb ) = -th( idxmb )/pi .* exp( -( s( idxmb ) + bump ).^2 ) .* -( s( idxmb ) + 1 ).^5 + exp( -( s( idxmb ) + bump ).^2 ) .* ( s( idxmb ) + 1 ).^4 .* ( 2 .* s( idxmb ).^2 + 12 .* s( idxmb ) + 5 );
% % out( idxob ) = -th( idxob )/pi .* -( s( idxob ) + 1 ).^5 - 5.* ( s( idxob ) + 1 ).^4;
% 
% out( idxpb ) = dp( th( idxpb ) );
% out( idxmb ) = dm( th( idxmb ) );
% out( idxob ) = do( th( idxob ) );
% 
% if( find( th == -pi ) )
%     out(  th == -pi  ) = 0;
% elseif( find( th == pi ) )
%     out(  th == pi  ) = 0;
% end

syms x
d = matlabFunction( diff( ( pi^2 - x^2 )/( 2 * pi ) * ( pi^2 - x^2 )^9 * -( 2 * atanh( x/pi ) + 1 )^10 * 2 * atanh( x/pi ) * ( 2 * atanh( x/pi ) - 1 )^4 * ( 2 * atanh( x/pi ) + 3 )^2, x ) );

out = d( th );
if( find( th == -pi ) )
    out(  th == -pi  ) = 0;
elseif( find( th == pi ) )
    out(  th == pi  ) = 0;
end