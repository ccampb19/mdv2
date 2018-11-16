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
            case {1,2}
                %Blue tape, fiber at 8 (or 10) o'clock
                NEWRHOS(ii,:) = sqrt((RHOS+0.45).^2+0.25.^2);
            case {3,6}
                %Green or red tape, fiber at 6 or 12 o'clock
                NEWRHOS(ii,:) = sqrt(RHOS.^2+.45^2);
            case {4,5}
                %Yellow or orange, fiber at 2 or 4 o'clock
                NEWRHOS(ii,:) = sqrt((RHOS-0.45).^2+0.25^2);
            otherwise
                error('Something Wrong...')
        end
end