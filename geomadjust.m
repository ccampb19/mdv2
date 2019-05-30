function NEWRHOS = geomadjust(USEDIODES,RHOS)
%MDLOAD loads files for multirho fits
%   RHOS contains nominal source-detector separations.
%   NEWRHOS are adjusted for ferrule geometry.
%   At this point, geometry is hard-coded in this function.

NEWRHOS = zeros(length(USEDIODES),length(RHOS));
for ii = 1:6
        switch ii
%             case {1}
                %B/W ;marking tape, center fiber: do nothing
            case {2,3}
                %Green or blue tape, fiber at 8 or 10 o'clock
                NEWRHOS(ii,:) = sqrt((RHOS+0.433).^2+0.25.^2);
            case {1,4}
                %Brown or yellow tape, fiber at 6 or 12 o'clock
                NEWRHOS(ii,:) = sqrt(RHOS.^2+.5^2);
            case {5,6}
                %Orange or red tape, fiber at 2 or 4 o'clock
                NEWRHOS(ii,:) = sqrt((RHOS-0.433).^2+0.25^2);
            otherwise
                error('Something Wrong...')
        end
end