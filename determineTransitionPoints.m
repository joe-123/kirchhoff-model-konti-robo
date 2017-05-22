function p = determineTransitionPoints(tcr, beta)
% Funktion gibt eine Matrix zurück.
% 1. Zeile steht für die transitionPoints (bzw. das entsprechende s)
% 2. Zeile beschreibt tube1 in den jeweiligen Segmenten
% 3. Zeile beschreibt tube2 ...
% --> 0=tube nicht vorhanden
% --> 1=tube gerade
% --> 2=tube gerkuemmt

    p = zeros(3, 0);
    % transition points
    for c = 1:length(beta)
        p1 = -beta(c) + tcr(c).tube.Ls;
        p2 = -beta(c) + tcr(c).tube.L;
        if (p1 > 0 && all(p(1,:)~=p1)) % p1 hinzufügen wenn nicht schon vorhanden
            p(1, end+1) = p1;
        end
        if (p2 > 0 && all(p(1,:)~=p2)) % p2 hinzufügen wenn nicht schon vorhanden
            p(1, end+1) = p2;
        end
    end
    
    % 0 hinzufügen wenn nicht schon vorhanden
    if all(p(1,:) ~= 0)
        p(1, end+1) = 0;
    end
    
    p = sort(p, 2);
    
    % status der tubes in einzelnen Segmenten --> 0=nicht vorhanden , 1=gerade , 2=gekrümmt
    for x = 1:length(p)
        if p(1, x) <= (-beta(2) + tcr(2).tube.L)
        
            % form von tube1
            if p(1, x) <= (-beta(1) + tcr(1).tube.Ls)
                p(2, x) = 1; % tube1 ist gerade
            else
                p(2, x) = 2; % tube1 is gekrümmt
            end
            
            % form von tube2
            if p(1, x) <= (-beta(2) + tcr(2).tube.Ls)
                p(3, x) = 1; % tube2 ist gerade
            else
                p(3, x) = 2; % tube2 is gekrümmt
            end
            
        else
            % form von tube1
            if p(1, x) <= (-beta(1) + tcr(1).tube.Ls)
                p(2, x) = 1; % tube1 ist gerade
            else
                p(2, x) = 2; % tube1 is gekrümmt
            end
            p(3, x) = 0; % tube 2 nicht vorhanden
        end
    end
    % Status der tubes muss in der matrix 1x nach links geschoben werden
    p(2:3, 1:end-1) = p(2:3, 2:end);
    p(2:3, end) = 0;
end