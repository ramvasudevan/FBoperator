function out = sampvf( th )

bump = 5;

out = zeros( size( th ) );
s = 2 * atanh( th/pi );

% idxpb = find( s > bump );
% idxmb = find( s < -bump );
% idxob = find( s <= bump & s >= -bump );
% 
% out( idxpb ) = ( pi^2 - th( idxpb ).^2 )/( 2 * pi ) .* exp( -( s( idxpb ) - bump ).^2 ) .* -( s( idxpb ) + 0.5 );
% out( idxmb ) = ( pi^2 - th( idxmb ).^2 )/( 2 * pi ) .* exp( -( s( idxmb ) + bump ).^2 ) .* -( s( idxmb ) + 0.5 );
% out( idxob ) = ( pi^2 - th( idxob ).^2 )/( 2 * pi ) .* -( s( idxob ) + 0.5 );
% 
% if( find( th == -pi ) )
%     out(  th == -pi  ) = 0;
% elseif( find( th == pi ) )
%     out(  th == pi  ) = 0;
% end

out = ( pi^2 - th.^2 )/( 2 * pi ) .* ( pi.^2 - th.^2 ).^9 .* -( s + 1 ).^10 .* s .* ( s - 1 ).^4 .* ( s + 3 ).^2;
if( find( th == -pi ) )
    out(  th == -pi  ) = 0;
elseif( find( th == pi ) )
    out(  th == pi  ) = 0;
end