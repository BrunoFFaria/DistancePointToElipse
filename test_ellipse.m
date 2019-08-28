clc
close all

phi = 0;
t=0:1/(2*pi):2*pi;
t1 = 290;
for u=80:10:360
    
    for ang=0:360
        % if t1 == 90
        %     t1 =89.9
        % end
        % phi = (2*t1)
        phi = u;
        semimajor_axis = 4;
        semiminor_axis = 2;
        
        x0 = 15;
        y0 = 15;
        t1 = ang;
        
        % calculo da elipse
        x_t1 = x0 + cosd(phi).*semimajor_axis.*cos(t) - sind(phi).*semiminor_axis.*sin(t);
        y_t1 = y0 + sind(phi).*semimajor_axis.*cos(t) + cosd(phi).*semiminor_axis.*sin(t);
        
        x_n11 = x0 + cosd(phi).*semimajor_axis.*cosd(t1) - sind(phi).*semiminor_axis.*sind(t1);
        y_n11 = y0 + sind(phi).*semimajor_axis.*cosd(t1) + cosd(phi).*semiminor_axis.*sind(t1);
        
        semimajor_axis = 12;
        semiminor_axis = 11;
        
        % calculo do ponto
        x_n12 = x0 + cosd(phi).*semimajor_axis.*cosd(t1) - sind(phi).*semiminor_axis.*sind(t1);
        y_n12 = y0 + sind(phi).*semimajor_axis.*cosd(t1) + cosd(phi).*semiminor_axis.*sind(t1);
        
        [semimajor_axis, semiminor_axis, x0, y0, phi,k] = ellipse_fit(x_t1, y_t1);% alteraçao dos eixos
        
        % obtenção dos focus da ellipse
        a0 = semimajor_axis;
        b0 = semiminor_axis;
        c = sqrt ( a0^2 - b0^2);
       
        F1.x = x0 - cos(phi)*c;
        F1.y = y0 - sin(phi)*c;
        
        F2.x = x0 + cos(phi)*c;
        F2.y = y0 + sin(phi)*c;
        
        % calculo dos declives
        m1 = (y_n12 - F1.y) / (x_n12 - F1.x);
        m2 = (y_n12 - F2.y) / (x_n12 - F2.x);
        
        c1 = F1.y - m1*F1.x;
        c2 = F2.y - m2*F2.x;
        
           
        % obtenção dos parametros da recta que intersecta os pontos da elipse
        B = sqrt(m2^2+1);
        A1 = sqrt(m1^2+1);
        A2 = -sqrt(m1^2+1);
        
        bi1 = (c1*B - c2*A1)/(B - A1);
        bi2 = (c1*B - c2*A2)/(B - A2);
        
        mi1 = ( m1*B - m2*A1)/(B - A1);
        mi2 = ( m1*B - m2*A2)/(B - A2);
   
        % escolha da bisectris
        
            % calculo da recta que passa pelos focos
            m_FC = (F2.y - F1.y)/(F2.x - F1.x);
            b_FC = F1.y - m_FC*F1.x;
            
            % calculo da intersecção das bisectrizes com a recta que passa
            % pelos focos
            
            % calculo da primeira bisectriz
            xb_i1 = (b_FC - bi1)/(mi1 - m_FC);
            yb_i1 = mi1*xb_i1 + bi1;
            
            % calculo da segunda bisectriz
            xb_i2 = (b_FC - bi2)/(mi2 - m_FC);
            yb_i2 = mi2*xb_i2 + bi2;
            
            F_11 = sqrt((F1.x-xb_i1)^2 + (F1.y-yb_i1)^2);
            F_12 = sqrt((F2.x-xb_i1)^2 + (F2.y-yb_i1)^2);
            F_21 = sqrt((F1.x-xb_i2)^2 + (F1.y-yb_i2)^2);
            F_22 = sqrt((F2.x-xb_i2)^2 + (F2.y-yb_i2)^2);
            F_1t = F_11 + F_12;
            F_2t = F_21 + F_22;
            
            if(F_1t <= F_2t)
                m = mi1;
                b = bi1;
            else
                m = mi2;
                b = bi2;
            end
            
        clc    
        
            
        % obtenção de y1 para x1 = 0;
        x1 = 1;
        y1 = m*x1 + b;
        
        % obtenção de y2 para x2 = x_n12
        x2 = 2;
        y2 = m*x2 + b;
        % obtenção dos t's
        y_intersecao = m*(-10:1:10)+b;
        % -> calculo dos parametros
        U = (x1 - x2);
        V = (y1 - y2);
        
        
        alpha = (k(1)*U^2 + 2*k(2)*U*V + k(3)*V^2);
        beta = ( 2 * k(1) * x1 * U + 2*k(2)*(x1*V + y1*U) + 2*k(3)*V*y1 + 2*k(4)*U + 2*k(5)*V);
        gamma = k(1) * x1^2 + 2*k(2) * x1*y1 + k(3)*y1^2 + 2*k(4)*x1 + 2* k(5)*y1 + k(6);
        
        % obtenção dos t's
        t1 = (-beta + sqrt(beta^2 - 4 * alpha * gamma))/(2*alpha);
        t2 = (-beta - sqrt(beta^2 - 4 * alpha * gamma))/(2*alpha);
        
        % obtenção dos pontos de intesecção
        xi1 = x1 + U*t1;
        yi1 = y1 + V*t1;
        
        xi2 = x1 + U*t2;
        yi2 = y1 + V*t2;
        u
        ang
        d1 = sqrt( (x_n12 - xi1)^2 + (y_n12 - yi1)^2)
        d2 = sqrt( (x_n12 - xi2)^2 + (y_n12 - yi2)^2)
        
        
        d = sqrt( (x_n12 - x_n11)^2 + (y_n12 - y_n11)^2)
        pause(1e-6)
        
        plot(x_t1,y_t1,'b',xi1,yi1,'r*',xi2,yi2,'g*',x_n12,y_n12,'*');
        hold on
        plot(-10:1:10,y_intersecao)
        plot([F2.x x_n12],[F2.y y_n12],'g',[F1.x x_n12],[F1.y y_n12],'g',[xb_i1 x_n12],[yb_i1 y_n12],'b',[xb_i2 x_n12],[yb_i2 y_n12],'k')
        %plot(0:1:20,m_FC*(0:1:20)+b_FC,'k')
        %plot([xb_i1 x_n12],[yb_i1 y_n12],'k',[xb_i2 x_n12],[yb_i2 y_n12],'g')
        %legend(num2str(ang+phi*180/pi))
        axis([0 40 0 40])
        hold off
    end
end
