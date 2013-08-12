function dy = vdpvf(t,y)

scaling = 1;
bumpx = 2.9;
bumpy = 2.9;
expbump = 200;
vdpscale = 1;

dy = zeros( 2, 1 );
% heaviside version
% if( abs( y( 1 ) ) > 1 ) 
%     dy( 1 ) = 0;
% else
%     dy( 1 ) = scaling * y( 2 );
% end
% 
% if( abs( y( 2 ) ) > 1 ) 
%     dy( 2 ) = 0;
% else
%     dy( 2 ) = -(scaling*y(1) + scaling*y(2)*((scaling*y(1))^2 - 1));
% end

% if( y(1) > pi )
%     y(1) = pi;
% elseif( y(1) < -pi )
%     y(1) = -pi;
% end
% 
% if( y(2) > pi )
%     y(2) = pi;
% elseif( y(2) < -pi )
%     y(2) = -pi;
% end

helper = [ scaling * y( 2 ); -(scaling*y(1) + scaling*vdpscale*y(2).*((scaling*y(1)).^2 - 1)) ];

% different exponential version
dy( 1 ) = ( pi^2 - y(1)^2 ) * ( pi^2 - y(2)^2 ) * helper( 1 );
dy( 2 ) = ( pi^2 - y(1)^2 ) * ( pi^2 - y(2)^2 ) * helper( 2 );

% exponential version
% dy(1) = exp(-1./(1-(y(1)/pi).^100)) .* ( helper( 1 ) + ( exp( -1 ) - exp(-1./(1-(y(2)/pi).^100) ) ) .* helper( 2 ) ) * exp( 1 );
% dy(2) = exp(-1./(1-(y(2)/pi).^100)) .* ( helper( 2 ) + ( exp( -1 ) - exp(-1./(1-(y(1)/pi).^100) ) ) .* helper( 1 ) ) * exp( 1 );

% dy(1) = exp(-1./(1-(y(1)/pi).^100)) .* exp(-1./(1-(y(2)/pi).^100)) .* ( helper( 1 ) ) * exp(1)^2;
% dy(2) = exp(-1./(1-(y(1)/pi).^100)) .* exp(-1./(1-(y(2)/pi).^100)) .* ( helper( 2 ) ) * exp(1)^2;


% dy(1) = scaling*y(2).*exp(-1./(1-(y(1)/pi).^100)) .* exp(-1./(1-(y(2)/pi).^100)) * (exp(1))^2;
% dy(2) = -(scaling*y(1) + scaling*y(2)*((scaling*y(1))^2 - 1))*exp(-1./(1-(y(1)/pi).^100)) .* exp(-1./(1-(y(2)/pi).^100)) * (exp(1))^2;