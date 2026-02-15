function [rhosKN,rhosSB,rhosOK,rhosUK] = getDike(rho1,rho2,y,thick)
%
% Initialisieren des Scripts.
na=length(y);
[rhosKN,rhosSB,rhosOK,rhosUK]=deal(zeros(size(y)));
anordnung = {'kleineNormale', 'Schlumberger', 'Oberkante', 'Unterkante'};

% Iteration über Konfigurationen.
for kk = 1:length(anordnung)
    cur_anordnung = anordnung{kk};
    
    % Festlegen der Elektrodenabstände
    switch cur_anordnung
        case 'kleineNormale'
            am=0.1;

        case {'Schlumberger', 'Oberkante', 'Unterkante'}
            am=0.4;
            mn=0.1;
    end
    
    % Berechnung der Potentiale an den Potentialsonden und des scheinbaren
    % spezifischen Widerstands.
    for ii=1:na
        switch cur_anordnung
            case 'kleineNormale'
                um=dike2(1,rho1,rho2,thick,0,y(ii)+am/2,0,y(ii)-am/2);
                rhosKN(ii)=um*2*pi./(1./am);

            case 'Schlumberger'
                uam=dike2(1,rho1,rho2,thick,0,y(ii)+(1)*(am+mn/2),0,y(ii)+(1)*mn/2);
                uan=dike2(1,rho1,rho2,thick,0,y(ii)+(1)*(am+mn/2),0,y(ii)+(1)*(-mn/2));
                ubm=dike2(1,rho1,rho2,thick,0,y(ii)-1*(am+mn/2),0,y(ii)-1*mn/2);
                ubn=dike2(1,rho1,rho2,thick,0,y(ii)-1*(am+mn/2),0,y(ii)-1*(-mn/2));
                % SB entspricht dem arithmetischen Mittel der OK und UK.
                rhosSB(ii)=0.5 * ((uam-uan)*2*pi./(1./am-1./(am+mn)) + ...
                                  (ubm-ubn)*2*pi./(1./am-1./(am+mn)));

            case 'Oberkante'
                um=dike2(1,rho1,rho2,thick,0,y(ii)+(1)*(am+mn/2),0,y(ii)+(1)*mn/2);
                un=dike2(1,rho1,rho2,thick,0,y(ii)+(1)*(am+mn/2),0,y(ii)+(1)*(-mn/2));
                rhosOK(ii)=(um-un)*2*pi./(1./am-1./(am+mn));

            case 'Unterkante'
                um=dike2(1,rho1,rho2,thick,0,y(ii)-1*(am+mn/2),0,y(ii)-1*mn/2);
                un=dike2(1,rho1,rho2,thick,0,y(ii)-1*(am+mn/2),0,y(ii)-1*(-mn/2));
                rhosUK(ii)=(um-un)*2*pi./(1./am-1./(am+mn));

        end
    end
end
