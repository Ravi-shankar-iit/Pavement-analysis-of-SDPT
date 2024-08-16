tic
syms x z zi et 
a=1;   b=1;                       % Dimension of plate
h=0.5;                            % Thickness of the plate
w=100;                            % local force vector for uniformly distributed load w
kw=0;         kp=0;               % Coefficients for Winkler and Pasternak foundation
load_case=input('Provide 1 for sinusoidal else any other number=')
%w=1*sin(pi*X/a)*sin(pi*Y/b);     % Sinusoidal load vector
Vload_node=[1];                   % Node numbers which are subjected to vehicle load
Vload=[00];                       % Transverse load at each node of line 9
a1=1;a2=1125;a12=0;               % Thermal expansion coefficients 
b1=0;b2=0.006;b12=0;              % Moisture expansion coefficients 

disp(' enter 1 to analyse for simply supported plate from all sides ')
disp(' enter 2 to analyse for clamped from all sides ') 
disp(' enter 3 to analyse for two edges simply supported two edges x0 and xa and other two free ')
disp(' enter 4 to analyse for  two edges simply supported two edges y0 and yb and  other two free ')
disp('enter 5 to analyse for simply supported from 2 sides (x0 and xa) agnd fixed from other two sides(y0 & yb)')
disp('enter 6 to analyse for simply supported from 2 sides (y0 and yb) and fixed from other two sides(x0 & xa)')
choicebc=1;                       % Selected boundary condition
m=14;   n=14;                       % Number of elements in x and y-directions(mesh size)
hx=a/m;     hy=b/n;               % Size of individual element
hxhalf=hx*0.5;
hyhalf=hy*0.5;
delT= 0*2*z/h;                         % Temperature variation function
%delM= 0.755462185*(z+h/2)/h;           % Moisture variation function
delM= 0*0.377731093*(z+h/2)/h;           % Moisture variation function
nl=1;                                  % Number of layers
h_layers = [1]*h;           % Height of Layers
Th=[0];                        % Angle of fibers in each layer
%% Material Properties 
LIBa=[21];                 % library number of each layer
h_center = sum(h_layers)/2;      % Calculate center height
Ha = zeros(2,nl);                % Initialize output matrix
% Calculate bottom and top surface for each layer
for i=1:nl
    z_bottom = h_center - sum(h_layers(1:i));
    z_top = h_center - sum(h_layers(1:i-1));
    Ha(1,i) = z_bottom;
    Ha(2,i) = z_top;             % Ha is defined top to bottom layer, with first row having bottom surface
end
 
DbarB=zeros(6);
DbarS=zeros(2);
N_thermal=zeros(6,1);
N_moisture=zeros(6,1);

for i=1:nl
    LIB=LIBa(i);
    h1=Ha(1,i); 
    h2=Ha(2,i);                  %h1 is bottom face, h2 is top
    hnl=h2-h1;
    %Material properties
    %         if LIB==1    
    %                 E11e= 30*10^9; E22= 7*10^9;  E33=0;  vx=0.15;   G12e=E11e/(2*(1+vx));  
    %                 G23=G12e;   G13=G23;    vx1=0.25;                       %PQC Layer 
    %          elseif LIB==2       
    %                 E11e= 7*10^9; E22= 7*10^9;  E33=0;  vx=0.15;   G12e=E11e/(2*(1+vx));  
    %                 G23=G12e;   G13=G23;    vx1=0.15;                       % DLC Layer 
    %         elseif LIB==3
    %                 E11e= 2*10^9; E22= 2*10^9;  E33=0;  vx=0.35;   G12e=E11e/(2*(1+vx));  
    %                 G23=G12e;   G13=G23;    vx1=0.35;                       %Granular Subbase
    %         elseif LIB==4
    %                 E11e= 5*10^7; E22= 5*10^7;  E33=0;  vx=0.35;   G12e=E11e/(2*(1+vx));  
    %                 G23=G12e;   G13=G23;    vx1=0.35;                       % Subgrade
    %         elseif LIB==5
    %                E11e=181*10^9; E22=10.3*10^9;  E33=0;  vx=0.28;   G12e= 7.17*10^9;  
    %                G13=G12e;   G23=2.39*10^9;    vx1=vx*E22/E11e;           % Hygrothermal
    %         elseif LIB==6
    %                E11e=66.2*10^9; E22=66.2*10^9;  E33=0;  vx=1/3;   G12e=E11e/(2*(1+vx));  
    %                G23=G12e;   G13=G23;    vx1=vx*E22/E11e;                 % metal(Titanium) 
    %         elseif LIB==7
    %                E11e=11.7*10^9; E22=11.7*10^9;  E33=0;  vx=1/3;   G12e=E11e/(2*(1+vx));  
    %                G23=G12e;   G13=G23;    vx1=vx*E22/E11e;                 % ceremic     
    %         elseif LIB==8
    %                E11e=181*10^9; E22=10.3*10^9;  E33=0;  vx=0.28;   G12e= 7.17*10^9;  
    %                G23=G12e;   G13=G23;    vx1=vx*E22/E11e;            
    %         elseif LIB==9    
    %             E11e=181*10^9; E22=10.3*10^9;  E33=0;  vx=0.28;   G12e= 7.17*10^9;  
    %             G13=G12e;   G23=2.39*10^9;    vx1=vx*E22/E11e; 
    %         elseif LIB==10
    %             E11e=25*10^9; E22=1*10^9;  E33=1*10^9;  vx=0.25;   G12e= 0.5*E22;  %Thermal
    %             G13=G12e;   G23=0.2*E22;    vx1=vx*E22/E11e;
    %         elseif LIB==11
    %             E11e=25*10^6; E22=1*10^6;  E33=1*10^9;  vx=0.25;   G12e= 0.5*E22;  
    %             G13=G12e;   G23=0.2*E22;    vx1=vx*E22/E11e;
    %         elseif LIB==12
    %             E11e=20.83; E22=10.94;  E33=10;  vx=0.44;   G12e= 6.1;  
    %             G13=3.71;   G23=6.19;    vx1=vx*E22/E11e;
    %         elseif LIB==13
    %             E11e=37.8; E22=13.1;  E33=13.1;  vx=0.25;   G12e= 8;  
    %             G13=8;   G23=5.3;    vx1=0.25;
    %         elseif LIB==14
    %             E11e=10^9; E22=10^9;  E33=10^9;  vx=0.3;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==15
    %             E11e=2.95*10^9; E22=2.95*10^9;  E33=10^9;  vx=0.35;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e; 
    %         elseif LIB==16
    %             E11e=3.4*10^9; E22=3.4*10^9;  E33=3.4*10^9;  vx=0.35;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==17
    %             E11e=3.025977861*10^9; E22=3.025977861*10^9;  E33=3.025977861*10^9;  vx=0.35;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==18
    %             E11e=23.12233278*10^9; E22=E11e;  E33=E22;  vx=0.3;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==19  % m_infinite
    %             E11e=2.959853509*10^9; E22=2.959853509*10^9;  E33=2.959853509*10^9;  vx=0.34;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==20  % m_infinite/2
    %             E11e=3.152067867*10^9; E22=3.152067867*10^9;  E33=3.152067867*10^9;  vx=0.34;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;
    %         elseif LIB==21  % m = 0
    %             E11e=3.4*10^9; E22=3.4*10^9;  E33=3.4*10^9;  vx=0.34;   G12e= E11e/(2*(1+vx));  
    %             G13=G12e;   G23=G13;    vx1=vx*E22/E11e;    
    %         end
    % 
    % %Constitutive material coefficients without elastic modulus        
    % C11e=1/(1-(vx*vx1));  C22=1/(1-(vx*vx1));  C33=0;  C44=1;    C55=(5/6);    C66=(5/6);
    % C12e=vx1/(1-(vx*vx1));   C21=vx/(1-(vx*vx1));   C13=0;   C23=0;    C31=0;    C32=0;
    % C14=0;  C15=0;  C16=0;  C24=0;  C25=0;  C26=0;  C34=0;  C35=0;  C36=0;
    % C41=0;  C42=0;  C43=0;  C45=0;  C46=0;
    % C51=0;  C52=0;  C53=0;  C54=0;  C56=0;
    % C61=0;  C62=0;  C63=0;  C64=0;  C65=0;
    
    Property_num = 1;                            % Prperty_num = 1,2,3,4    (Em = 3.113391309   ( 0.4, 0.5, 0.6, 0.7)) 
                                                 % Prperty_num = 5,6,7,8    (Em = 3.245923216   ( 0.4, 0.5, 0.6, 0.7)) 
    Q_effective = Effective_prop (Property_num); % Prperty_num = 9,10,11,12 (Em = 3.4           ( 0.4, 0.5, 0.6, 0.7))   

     %Final constitutive coefficients for no angle ply
    Q11=Q_effective(1,1);   Q12=Q_effective(1,2);   Q13=Q_effective(1,3);    Q21=Q_effective(2,1);    Q22=Q_effective(2,2);    Q23=Q_effective(2,3);
    Q31=Q_effective(3,1);   Q32=Q_effective(3,2);   Q33=Q_effective(3,3);    Q44=Q_effective(4,4);    Q55=Q_effective(5,5);    Q66=Q_effective(6,6);
    Q14=0;  Q15=0;  Q16=0;  Q24=0;  Q25=0;  Q26=0;  Q34=0;  Q35=0;  Q36=0;
    Q41=0;  Q42=0;  Q43=0;  Q45=0;  Q46=0;
    Q51=0;  Q52=0;  Q53=0;  Q54=0;  Q56=0;
    Q61=0;  Q62=0;  Q63=0;  Q64=0;  Q65=0;
    
    % Relations for fibre angle profile
    Q11r=Q11*(cos(Th(i)))^4+2*(Q12+2*Q44)*(sin(Th(i)))^2*(cos(Th(i)))^2+Q22*(sin(Th(i)))^4;
    Q12r=(Q11+Q22-4*Q44)*(sin(Th(i)))^2*(cos(Th(i)))^2+Q12*((sin(Th(i)))^4+(cos(Th(i)))^4);
    Q22r=Q11*(sin(Th(i)))^4+2*(Q12+2*Q44)*(sin(Th(i)))^2*(cos(Th(i)))^2+Q22*(cos(Th(i)))^4;
    Q14r=(Q11-Q12-2*Q44)*(sin(Th(i)))*(cos(Th(i)))^3+(Q12-Q22+2*Q44)*(sin(Th(i)))^3*cos(Th(i));
    Q24r=(Q11-Q12-2*Q44)*(sin(Th(i)))^3*(cos(Th(i)))+(Q12-Q22+2*Q44)*(sin(Th(i)))*(cos(Th(i)))^3;
    Q44r=(Q11+Q22-2*Q12-2*Q44)*(sin(Th(i)))^2*(cos(Th(i)))^2+Q44*((sin(Th(i)))^4+(cos(Th(i)))^4);
    Q55r=Q55*(cos(Th(i)))^2+Q66*(sin(Th(i))^2);
    Q56r=(Q66-Q55)*cos(Th(i))*sin(Th(i));
    Q66r=Q55*(sin(Th(i)))^2+Q66*(cos(Th(i)))^2;
   
    ax=a1*cos(Th(i))^2+a2*sin(Th(i))^2; ay=a1*sin(Th(i))^2+a2*cos(Th(i))^2;     axy=2*(a1-a2)*sin(Th(i))*cos(Th(i));
    bx=b1*cos(Th(i))^2+b2*sin(Th(i))^2; by=b1*sin(Th(i))^2+b2*cos(Th(i))^2;     bxy=2*(b1-b2)*sin(Th(i))*cos(Th(i));
    a_thermal(:,i) =[ax;ay;axy];                                            % Thermal coefficients for each layer
    a_moisture(:,i)=[bx;by;bxy];                                            % Moisture coefficients for each layer
    Cg(:,:,i)=[Q11r Q12r Q14r;Q12r Q22r Q24r;Q14r Q24r Q44r];               % Bending constitutive terms
    CSg(:,:,i)=[Q55r  Q56r; Q56r Q66r];                                     % Transverse shear constitutive terms          
    CB=Cg(:,:,i);    CS=CSg(:,:,i);
    Gb=[1 0 0 z 0 0; 0 1 0 0 z 0; 0 0 1 0 0 z];
    
    DbarB=DbarB+ vpa(int(Gb'*CB*Gb,z,h1,h2));
    DbarS=DbarS+ vpa(int(CS,z,h1,h2));
    NT=(Gb'*CB*[ax;ay;axy]*delT);                                          %Thermal coefficients for each layer without z-integral
    NM=(Gb'*CB*[bx;by;bxy]*delM);                                          %Moisture coefficients for each layer without z-integral
    N_thermal=N_thermal+int(NT,z,h1,h2);
    N_moisture=N_moisture+int(NM,z,h1,h2);
end
%% stiffness matrix calculation
phi1=0.25*(1-zi)*(1+et)*(-zi+et-1);
phi2=0.5*(1-(zi^2))*(1+et);
phi3=0.25*(1+zi)*(1+et)*(et+zi-1);
phi4=0.5*(1+zi)*(1-(et^2));
phi5=0.25*(1+zi)*(1-et)*(zi-et-1);
phi6=0.5*(1-(zi^2))*(1-et);
phi7=0.25*(1-zi)*(1-et)*(-zi-et-1);
phi8=0.5*(1-zi)*(1-(et^2));
phi=[phi1 phi2 phi3 phi4 phi5 phi6 phi7 phi8];
no=size(phi',1);
%             
%       7-----6------5
%       |            |
%       |            |
%       8      ------4---->
%       |      |
%              |     |
%       |      |     |
%       1------2------3
%              |
%              \/

Ono=zeros(1,8);
%% Jacobian and dphi/dx,dphi/dy
x=[-hxhalf 0 hxhalf hxhalf hxhalf 0 -hxhalf -hxhalf];
y=[hyhalf hyhalf hyhalf 0 -hyhalf -hyhalf -hyhalf 0];
j=[diff(phi,zi);diff(phi,et)]*[x' y'];
detJ=det(j); %j(1,1)*j(2,2)-j(1,2)*j(2,1);
j11=j(2,2)/detJ;
j12=-j(1,2)/detJ;
j21=-j(2,1)/detJ;
j22=j(1,1)/detJ;
dpx=vpa(j11*diff(phi,zi)+j12*diff(phi,et));
dpy=vpa(j21*diff(phi,zi)+j22*diff(phi,et));
dp2x=vpa(j11*diff(dpx,zi)+j12*diff(dpx,et));
dp2y=vpa(j21*diff(dpy,zi)+j22*diff(dpy,et));

BbarB=[dpx,Ono,Ono,Ono,Ono;
       Ono,dpy,Ono,Ono,Ono;
       dpy,dpx,Ono,Ono,Ono;
       Ono,Ono,Ono,dpx,Ono;
       Ono,Ono,Ono,Ono,dpy;
       Ono,Ono,Ono,dpy,dpx];

BbarS=[Ono,Ono,dpx,phi,Ono;
       Ono,Ono,dpy,Ono,phi];
%% Formulation of element stiffness matrix/load vector
NOE=m*n;                %total number of elements
TONN=3*m*n+2*n+2*m+1;   %total number of nodes
kg=zeros(TONN*5);       %TOTAL STIFFNESS MATRIX
fg=zeros(5*TONN,1);     %Total load vector
ft_g=zeros(5*TONN,1);   %Total thermal load vector
fm_g=zeros(5*TONN,1);   %Total moisture load vector

i=[1:5]';
%av,bv,c,d,e,f,g,hv ..  are Dof number corresponsing to node number n1,n2,n3,
% n4,n5,n6,n7 and n8, respectively 
%These are for element
  av(i)=8*i-7;
  bv(i)=8*i-6; 
  c(i)=8*i-5;
  d(i)=8*i-4;
  e(i)=8*i-3;
  f(i)=8*i-2;
  g(i)=8*i-1;
  hv(i)=8*i;
  E=[av' bv' c' d' e' f' g' hv'];

% Defintion of element node connectivity
% element connectivitry matrix for saving  element number with node number
elementcon=zeros(NOE,9); %first column element number
ke=zeros(40,40,NOE);
fl=zeros(8,NOE);    %Contribution of Mechanical load to w DODFs only
ft=zeros(40,NOE);   %Contribution of Thermal to all DODFs 
fm=zeros(40,NOE);   %Contribution of moisture to all DODFs

parfor ne=1:NOE
    tic
%for ne=1:NOE
    q=floor(ne/m);              %element row number q+1
    r=rem(ne,m);                %element number in (q+1)th row
    if r==0
        n1=2*(q-1)+1;
        n2=n1+1;
        n3=n1+2 ;
        node1=((n1-1)*(2*m+1)/2)+((n1-1)*(m+1)*0.5)+(2*m-1);
        node2=node1+1;
        node3=node1+2;
        node7=((n3-1)*(2*m+1)*0.5)+((n3-1)*(m+1)*0.5)+(2*m-1);
        node6=node7+1;
        node5=node7+2;
        node8=(n2*(2*m+1)*0.5)+((n2-2)*(m+1)*0.5)+m;
        node4=node8+1;
        N=[node1 node2 node3 node4 node5 node6  node7 node8]';
        X=(a-hxhalf)+(zi*hxhalf);       Y=(q-1)*hy+hyhalf-(et*hyhalf);      %Point on the element for load definition
    else
        n1=2*q+1;
        n2=n1+1;
        n3=n1+2;
        node1=((n1-1)*(2*m+1)/2)+((n1-1)*(m+1)*0.5)+(2*r-1);
        node2=node1+1;
        node3=node1+2;
        node7=((n3-1)*(2*m+1)*0.5)+((n3-1)*(m+1)*0.5)+(2*r-1);
        node6=node7+1;
        node5=node7+2;
        node8=(n2*(2*m+1)*0.5)+((n2-2)*(m+1)*0.5)+r;
        node4=node8+1;
        N=[node1 node2 node3 node4 node5 node6  node7 node8]';
        X=(r-1)*hx+hxhalf+(zi*hxhalf);      Y=q*hy+hyhalf-(et*hyhalf);      %Point on the element for load definition
    end
    elementcon(ne,:)=[ne, N'];
    coordinate1(ne,1)=X;   coordinate2(ne,1)=Y;
    if load_case==1
        fun_sinusoidal=sin(pi*X/a)*sin(pi*Y/b);
    else
        fun_sinusoidal=1;
    end
    flv=w*fun_sinusoidal*phi';                                              %Case for UDL or sinusoidal load over the plate
    
 %Three point gauss integration for bending stiffness and load vector
    vb=[-sqrt(3/5) 0 sqrt(3/5) ];                                           %Gauss points
    wb=[0.5555555555 0.8888888889 0.5555555555];                            %Weight functions
    
    ab=0;   bb=0;                                                           %Initialization
    for i=1:3
        for j=1:3
            ab=vb(i);
            bb=vb(j);
            BbarB1=subs(subs(BbarB,zi,ab),et,bb);
            Kb=BbarB1'*DbarB*BbarB1;        
            keB=vpa(Kb);  
            F_thermal=BbarB1'*N_thermal*fun_sinusoidal;
            F_moisture=BbarB1'*N_moisture*fun_sinusoidal;
            K_winkler=kw*phi'*phi;      K_paster=kp*(dpx'*dpx+dpy'*dpy);    %Contribution terms for elastic foundation part
            K_elastic=subs(subs((K_winkler+K_paster),zi,ab),et,bb);         %Contribution of elastic foundation at Gauss points
            ke(:,:,ne)=ke(:,:,ne)+ vpa(keB*wb(i)*wb(j)*subs(subs(detJ,zi, ab),et,bb));      %Stiffness term after integral
            k_el=vpa(K_elastic*wb(i)*wb(j)*subs(subs(detJ,zi, ab),et,bb));  %Stiffness due to elastic foundation after integral (size 8x8)
            %ke(17:24,17:24)=ke(17:24,17:24)+k_el; Uncomment if winkler's is used.                       %Total stiffness contribution
            fl(:,ne)=fl(:,ne)+ vpa(subs(vpa(subs(flv,zi, ab)),et,bb)*wb(i)*wb(j)*subs(subs(detJ,zi, ab),et,bb)); %only for w direction (size 8x1)
            ft(:,ne)=ft(:,ne)+ vpa(wb(i)*wb(j)*subs(subs(F_thermal*detJ,zi, ab),et,bb)); %Thermal load contribution included
            fm(:,ne)=fm(:,ne)+ vpa(wb(i)*wb(j)*subs(subs(F_moisture*detJ,zi, ab),et,bb)); %Moisture load contribution included
        end
    end
  
    % Two point gauss integration for shear stiffness
    vs=[-sqrt(1/3) sqrt(1/3) ];     %Gauss points
    ws=[1 1];                       %Weight functions
    
    for i=1:2
    %parfor i=1:2
        for j=1:2
            as=vs(i);
            bs=vs(j);
            BbarS1=subs(subs(BbarS,zi,as),et,bs);
            Ks=BbarS1'*DbarS*BbarS1;
            keS=vpa(Ks);
            ke(:,:,ne)=ke(:,:,ne)+ vpa(keS*ws(i)*ws(j)*subs(subs(detJ,zi, as),et,bs));
        end
    end
    toc
end
%Order of Ke is u01..u08,v01...v08,---
%fl define load vector for w displacement only at 8 nodes
%Order of ft is u01..u08,v01...v08,---
%% Assembly of global stiffness matrix and load vector    
%Fixing of 8x8 matrix with 5 dof at each node in global kg matrix
  for ne = 1:NOE
    for i=1:8
        for j=1:8
            aa1=5*elementcon(ne,i+1)-4;     %for node number N(i) 1st Dof
            aa2=5*elementcon(ne,i+1);       %for node number N(i) last Dof
            bb1=5*elementcon(ne,j+1)-4;
            bb2=5*elementcon(ne,j+1);
            ee1=E(:,i);
            ee2=E(:,j);
            kg(aa1:aa2,bb1:bb2)= kg(aa1:aa2,bb1:bb2)+ke(ee1,ee2,ne);
        end
            ft_g(aa1:aa2,1)=ft_g(aa1:aa2,1)+ft(ee1,ne);
            fm_g(aa1:aa2,1)=fm_g(aa1:aa2,1)+fm(ee1,ne);
    end

    %Order of kg is with respect to nodes. so DOFs of each node are
    %numbered completely one by one  5+5+5+....total nodes
   
    %taking positions vertical displacement (w)DoF
        ww=5*elementcon(ne,2:9)'-2; %DOUBT
        fg(ww,1)=fg(ww,1)+fl(:,ne); %loading vector corresponding to the transverse direction
 end
%% Load contribution to specific nodes
Vload_dof=5*Vload_node-2*ones(size(Vload_node,1),1);
fg(Vload_dof)=fg(Vload_dof)+Vload;

%% Boundary condition
R=Bound_cond(m,choicebc);
%% Pliminary Results
adof=[1:(5*TONN)]';
adof(R,:)=[];
kga=kg(adof,adof);
fga=fg(adof,1)+ft_g(adof,1)+fm_g(adof,1);
uu=kga\fga;
uf=zeros(5*TONN,1);
uf(adof,1)=uf(adof,1)+uu;
noi=(floor(TONN/2)+1)*5-2;
unoi=uf(noi);
% non dimensional displacement
%dlnoi=unoi*E2F*h^3*10/(w*a^4)
%% Post processing of Stresses
% Node number sharing the mid node
% elementno = m * 0.5 * (n - 1);
element_range = m*((m/2)-1)+1:m*(m/2); % Define the range of element numbers
%element_range = 1:m;
num_elements = length(element_range);

% Initialize the 3D matrix to store stresses
%sz1 = zeros(3, 9, num_elements);
stresses_3d_xx = zeros(9, nl, num_elements, 'sym');
stresses_3d_yy = zeros(9, nl, num_elements, 'sym');
stresses_3d_xy = zeros(9, nl, num_elements, 'sym');
stresses_3d_xz = zeros(4, nl, num_elements, 'sym');
stresses_3d_yz = zeros(4, nl, num_elements, 'sym');

% Loop through the element numbers
for elem_idx = 1:num_elements
    elementno = element_range(elem_idx);
    nodenumber = elementcon(elementno, 2:9); % node numbers in order for given element
    X = coordinate1(elementno, 1);
    Y = coordinate2(elementno, 1);

    if load_case == 1
        fun_sinusoidal = sin(pi * X / a) * sin(pi * Y / b);
    else
        fun_sinusoidal = 1;
    end

    % Specific displacement parameters corresponding to all nodes of given element
    u0ele = uf((5 * nodenumber) - 4 * ones(1, 8))';
    theta1ele = uf((5 * nodenumber) - ones(1, 8))';
    v0ele = uf((5 * nodenumber) - 3 * ones(1, 8))';
    theta2ele = uf((5 * nodenumber))';
    woele = uf((5 * nodenumber) - 2 * ones(1, 8))';

    % Calculation of strains
    ep1 = j11 * diff(phi, zi) + j12 * diff(phi, et); % d/dx(phi)
    ep2 = j21 * diff(phi, zi) + j22 * diff(phi, et); % d/dy(phi)
    ep3 = phi;

    Exx = (ep1 * u0ele') + z * (ep1 * theta1ele');
    Eyy = (ep2 * v0ele') + z * (ep2 * theta2ele');
    Exy = (ep2 * u0ele' + ep1 * v0ele') + z * (ep2 * theta1ele' + ep1 * theta2ele');
    Exz=  (ep1 * woele')+ (ep3 * theta1ele') ;
    Eyz=  (ep2 * woele')+ (ep3 * theta2ele') ;
    inplane_strainZ    = [Exx; Eyy; Exy]; % function of zi, et and z
    transverse_strainZ = [Exz;Eyz];       % function of zi, et and z
    sigxx = zeros(9, nl, 'sym');
    sigxy = zeros(9, nl, 'sym');
    sigyy = zeros(9, nl, 'sym');
    sigxz = zeros(4, nl, 'sym');
    sigyz = zeros(4, nl, 'sym');
    
    for NL = 1:nl
        Strain_Thermal = a_thermal(:, NL) * delT * fun_sinusoidal;
        Strain_Moisture = a_moisture(:, NL) * delM * fun_sinusoidal;
        inplane_stressZ = Cg(:, :, NL) * (inplane_strainZ - Strain_Thermal-Strain_Moisture);
        transverse_stressZ=CSg(:,:,NL)*transverse_strainZ;                                        % multiplay it with 2 
        vb = [-sqrt(3/5) 0 sqrt(3/5)];

        for i = 1:3
            for j = 1:3
                ZI = vb(j);
                ET = vb(i);
                p = 3 * (i - 1) + j;
                inplane_stress_gauss(:, p) = subs(subs(inplane_stressZ, zi, ZI), et, ET);
            end
        end
        vs=[-sqrt(1/3) sqrt(1/3) ];     %Gauss points
        for i=1:2
            for j=1:2
                ZI=vs(j);   ET=vs(i);
                p=2*(i-1)+j;
                transverse_stress_gauss(:,p) = (6/5)*subs(subs(transverse_stressZ,zi,ZI),et,ET);
            end
        end
        % Accumulate stresses for each layer
        sigxx(:, NL) = inplane_stress_gauss(1, :)';
        sigyy(:, NL) = inplane_stress_gauss(2, :)';
        sigxy(:, NL) = inplane_stress_gauss(3, :)';
        sigxz(:, NL) = transverse_stress_gauss(1, :)';
        sigyz(:, NL) = transverse_stress_gauss(2, :)';
    end
    
    % Store stresses for each element
    stresses_3d_xx(:, :, elem_idx) = sigxx;
    stresses_3d_yy(:, :, elem_idx) = sigyy;
    stresses_3d_xy(:, :, elem_idx) = sigxy;
    stresses_3d_xz(:, :, elem_idx) = sigxz;
    stresses_3d_yz(:, :, elem_idx) = sigyz;
end
toc 
%% graph 
stresses_num_xx = stresses_3d_xx(:, :,m/2);
stresses_num_yy = stresses_3d_yy(:, :,m/2);
stresses_num_xy = stresses_3d_xy(:, :,m/2);
stresses_num_xz = stresses_3d_xz(:, :,1);
stresses_num_yz = stresses_3d_yz(:, :,1);
sigxxn = zeros(11,9,nl);
sigyyn = zeros(11,9,nl);
sigxyn = zeros(11,9,nl);
sigxzn = zeros(11,4,nl);
sigyzn = zeros(11,4,nl);

for layer = 1:nl
    h1 = Ha(1, layer);
    h2 = Ha(2, layer);
    increment = (h2 - h1) / 10;
    k = 1;
    
    for i = h1:increment:h2

              sigxxn(k, :, layer) = vpa(1*subs(stresses_num_xx(:, layer,:), z, i));
              sigyyn(k, :, layer) = vpa(1*subs(stresses_num_yy(:, layer,:), z, i));
              sigxyn(k, :, layer) = vpa(1*subs(stresses_num_xy(:, layer,:), z, i));
              sigxzn(k, :, layer) = vpa(1*subs(stresses_num_xz(:, layer,:), z, i));
              sigyzn(k, :, layer) = vpa(1*subs(stresses_num_yz(:, layer,:), z, i));
              % sigxxn(k, :, layer) = vpa(subs(stresses_3d_xx(:, :,2), z, i));
              % sigyyn(k, :, layer) = vpa(subs(stresses_3d_yy(:, :, layer, elem_idx)', z, i));
              % sigxyn(k, :, layer) = vpa(subs(stresses_3d_xy(:, :, layer, elem_idx)', z, i));
             
        depth(k, layer) = i;
        k = k + 1;
    end
%% RESULTS FOR PAPERS
sigma_xx=sigxxn(11,9,1); % At (a/2,b/2,h/2)
sigma_xz=sigxzn(6,2,1);

    figure (1)
    if((layer-1)~=0)
        %line 1
        disp("Hello");
        x = [sigxxn(1, 9, layer-1); sigxxn(end, 9, layer)];
        y = [depth(1, layer-1)/h; depth(end, layer)/h];
        plot(x, y, "LineStyle","-.","Color","k","LineWidth",0.7);
        legend("hide")
        hold on
    end
    plot(sigxxn(:,9,layer), depth(:, layer)/h, "LineStyle","-","Color","k","LineWidth",1,"Marker","*");
        xlabel("\sigma_x(N/mm^2)","FontSize",14);
        ylabel("z*","FontSize",14);
        legend("hide")
        hold on 
        figure (2)
        if((layer-1)~=0)
        %line 1
        disp("Hello");
        x = [sigyyn(1, 9, layer-1); sigyyn(end, 9, layer)];
        y = [depth(1, layer-1)/h; depth(end, layer)/h];
        plot(x, y,"LineStyle","-.","Color","k","LineWidth",0.7);
        legend("hide")
        hold on
        end
    plot(sigyyn(:,9,layer), depth(:, layer)/h, "LineStyle","-","Color","k","LineWidth",1,"Marker","*");
        xlabel("{\sigma}_{y}(N/mm^2)","FontSize",14);
        ylabel("z*","FontSize",14);
        legend("hide")
        hold on 
        figure (3)
       if((layer-1)~=0)
        %line 1
        disp("Hello");
        x = [sigxyn(1, 9, layer-1); sigxyn(end, 9, layer)];
        y = [depth(1, layer-1)/h; depth(end, layer)/h];
        plot(x, y,"LineStyle","-.","Color","k","LineWidth",0.7);
        legend("hide")
        hold on
    end
    plot(sigxyn(:,9,layer), depth(:, layer)/h, "LineStyle","-","Color","k","LineWidth",1,"Marker","*");
        xlabel("{\sigma}_{xy}(N/mm^2)","FontSize",14);
        ylabel("z*","FontSize",14);
        legend("hide")
        hold on 
            figure (4)
        if((layer-1)~=0)
        %line 1
        disp("Hello");
        x = [sigxzn(1, 4, layer-1); sigxzn(end, 4, layer)];
        y = [depth(1, layer-1)/h; depth(end, layer)/h];
        plot(x, y,"LineStyle","-.","Color","k","LineWidth",0.7);
        legend("hide")
        hold on
        end
    plot(sigxzn(:,4,layer), depth(:, layer)/h, "LineStyle","-","Color","k","LineWidth",1,"Marker","*");
        xlabel("{\sigma}_{xz}(N/mm^2)","FontSize",14);
        ylabel("z*","FontSize",14);
        legend("hide")
        hold on 
            figure (5)
        if((layer-1)~=0)
        %line 1
        disp("Hello");
        x = [sigyzn(1, 2, layer-1); sigyzn(end, 2, layer)];
        y = [depth(1, layer-1)/h; depth(end, layer)/h];
        plot(x, y,"LineStyle","-.","Color","k","LineWidth",0.7);
        legend("hide")
        hold on
    end
    plot(sigyzn(:,2,layer), depth(:, layer)/h, "LineStyle","-","Color","k","LineWidth",1,"Marker","*");
        xlabel("{\sigma}_{yz}(N/mm^2)","FontSize",14);
        ylabel("z*","FontSize",14);
        legend("hide")
        hold on 
end
% %% Excel Sheet data 
% % Define file name
% excelFileName = 'resin .xlsx';
% 
% % Initialize matrices to store data
% data_sigxx = zeros(11, nl);
% data_sigyy = zeros(11, nl);
% data_sigxy = zeros(11, nl);
% data_sigxz = zeros(11, nl);
% data_sigyz = zeros(11, nl);
% 
% for layer = 1:nl
%     % Extract data for current layer
%     data_sigxx(:, layer) = sigxxn(:,9,layer);
%     data_sigyy(:, layer) = sigyyn(:,9,layer);
%     data_sigxy(:, layer) = sigxyn(:,9,layer);
%     data_sigxz(:, layer) = sigxzn(:,4,layer);
%     data_sigyz(:, layer) = sigyzn(:,2,layer);
% end
% 
% % Write data to Excel file
% xlswrite(excelFileName, data_sigxx, 'sigxx');
% xlswrite(excelFileName, data_sigyy, 'sigyy');
% xlswrite(excelFileName, data_sigxy, 'sigxy');
% xlswrite(excelFileName, data_sigxz, 'sigxz');
% xlswrite(excelFileName, data_sigyz, 'sigyz');
%% along x axis 
layer_x=1;
%stresses_num_xx_liner = subs(stresses_3d_xx(7:9, layer_x,1:4),z,h/2);

stresses_num_xx_liner_squeeze = squeeze(subs(stresses_3d_xx(7:9, layer_x,1:m),z,h/2));
stresses_num_xx_liner_reshape = stresses_num_xx_liner_squeeze(:);
stresses_num_yy_liner_squeeze = squeeze(subs(stresses_3d_yy(7:9, layer_x,1:m),z,h/2));
stresses_num_yy_liner_reshape = stresses_num_yy_liner_squeeze(:);
stresses_num_xy_liner_squeeze = squeeze(subs(stresses_3d_xy(7:9, layer_x,1:m),z,h/2));

stresses_num_xy_liner_reshape = stresses_num_xy_liner_squeeze(:);
%global_x_cooridate = 0:a/(3*m-1):a;
%% Globle x coordinate 
% Define global coordinate system parameters
vb1=[-sqrt(3/5) 0 sqrt(3/5) ]; %(for parfor execution)
global_length = num_elements * hx;  % mm

% Initialize an array to store the global Gauss points
global_gauss_points = zeros(1, num_elements * numel(vb1));

% Index to keep track of the current position in the global_gauss_points array
global_gauss_index = 1;

% Calculate and store global Gauss points for all elements
for element_index = 0:(num_elements - 1)
    for local_x = vb1
        global_x = element_index * hx + (hx / 2) * (1 + local_x);
        global_gauss_points(global_gauss_index) = global_x;
        global_gauss_index = global_gauss_index + 1;
    end
end

% Now you have the global Gauss points stored in the 'global_gauss_points' array
%% 
% stresses_num_xx_liner_reshape = reshape(stresses_num_xx_liner_squeeze, size(stresses_num_xx_liner_squeeze, 1), []);
         midline_delfection = 1.5*m^2+m+1:1.5*m^2+3*m+1;
         X_axis = 0:a/(length(midline_delfection)-1):a;
         figure (6)
         plot ( X_axis, uf(5*midline_delfection -2), "LineStyle","--","Color","g","LineWidth",1.0,"Marker","*");
             xlabel(" Length(m) ");
             ylabel("Deflection (mm)");
             hold on
        figure (7)
             plot ( global_gauss_points, stresses_num_xx_liner_reshape,"LineStyle","-","Color","g","LineWidth",1.0,"Marker","*")
             xlabel(" Length(m) ");
             ylabel("{\sigma}_{xx}(N/m^2)","FontSize",14);
             hold on
        figure (8)
             plot ( global_gauss_points, stresses_num_yy_liner_reshape,"LineStyle","-","Color","g","LineWidth",1.0,"Marker","*")
             xlabel(" Length(m)");
             ylabel("{\sigma}_{yy}(N/m^2)","FontSize",14);
             hold on

        figure (9)
             plot ( global_gauss_points, stresses_num_xy_liner_reshape,"LineStyle","-","Color","g","LineWidth",1.0,"Marker","*")
             xlabel(" Length(m) ");
             ylabel("{\sigma}_{xy}(N/m^2)","FontSize",14);
             hold on

%% dta excel 
% % Define the arrays to store data points
% deflection_data = -1*uf(5*midline_delfection -2);
% stress_xx_data = stresses_num_xx_liner_reshape;
% stress_yy_data = stresses_num_yy_liner_reshape;
% stress_xy_data = stresses_num_xy_liner_reshape;
% 
% % Check dimensions
% if numel(X_axis) ~= numel(deflection_data) || ...
%    numel(global_gauss_points) ~= numel(stress_xx_data) || ...
%    numel(global_gauss_points) ~= numel(stress_yy_data) || ...
%    numel(global_gauss_points) ~= numel(stress_xy_data)
%     error('Dimensions of arrays are not consistent.');
% end
% 
% % Plotting the figures
% figure(6)
% plot(X_axis, deflection_data, "LineStyle","--","Color","b","LineWidth",1.0,"Marker","square");
% xlabel("Length (x)");
% ylabel("Deflection");
% hold on
% 
% figure(7)
% plot(global_gauss_points, stress_xx_data, "LineStyle","-","Color","b","LineWidth",1.0,"Marker","square")
% xlabel("Length (x)");
% ylabel("sigma xx");
% hold on
% 
% figure(8)
% plot(global_gauss_points, stress_yy_data, "LineStyle","-","Color","b","LineWidth",1.0,"Marker","square")
% xlabel("Length (x)");
% ylabel("sigma yy");
% hold on
% 
% figure(9)
% plot(global_gauss_points, stress_xy_data, "LineStyle","-","Color","b","LineWidth",1.0,"Marker","square")
% xlabel("Length (x)");
% ylabel("sigma xy");
% hold on
% 
% % Write data to Excel
% data_matrix = [X_axis(:), deflection_data(:), stress_xx_data(:), stress_yy_data(:), stress_xy_data(:)];
% header = {'Length (x)', 'Deflection', 'sigma xx', 'sigma yy', 'sigma xy'};
% writematrix(header.', 'output.xlsx', 'Sheet', 'Data', 'Range', 'A1');
% writematrix(data_matrix, 'output.xlsx', 'Sheet', 'Data', 'Range', 'A2');

%% ro
% sigxxn = zeros(11,9,nl);
% sigyyn = zeros(11,9,nl);
% sigxyn = zeros(11,9,nl);
% 
% for layer=1:nl
%     h1=Ha(1,layer); h2=Ha(2,layer);
%     increment=(h2-h1)/10;
%     k=1;
%     for i = h1:increment:h2
%         sigxxn(k,:,layer)= vpa(subs(sigxx(:,:,layer)',z,i));
%         sigyyn(k,:,layer)= vpa(subs(sigyy(:,:,layer)',z,i));
%         sigxyn(k,:,layer)= vpa(subs(sigxy(:,:,layer)',z,i));
%         depth(k)=i;
%         k= k+1;
%     end
%     figure (1)
%     plot(depth/h, 10^-6*sigxxn(:,9,layer));
%         xlabel("Depth/h");
%         ylabel("sigxx");
%    hold on
%     figure (2)
%     plot(depth/h,10^-6*sigxyn(:,9,layer));
%             xlabel("Depth");
%             ylabel("sigxy");
%     hold on
%     figure (3)
%     plot(depth/h,10^-6*sigyyn(:,9,layer));
%             xlabel("Depth");
%             ylabel("sigyy");
%     hold on
% 
%     figure (4)
%     plot(10^-6*sigxxn(:,9,layer),depth/h);
%             xlabel("sigyy");
%             ylabel("Depth");
%     hold on   
% 
%     midline = [29:37];
%     X_axis = 0:a/(length(midline)-1):a;
% 
%     figure (5)
%     plot ( X_axis, -1*uf(5*midline -2));
%             xlabel("X axis ");
%             ylabel("Deflection");
% end
 %%
% for layer = 1:nl
%     h1 = Ha(1, layer);
%     h2 = Ha(2, layer);
%     increment = (h2 - h1) / 10;
%     k = 1;
% 
%     for i = h1:increment:h2
% 
%               sigxxn(k, :, layer) = vpa(subs(stresses_num_xx(:, layer,:), z, i));
%               sigyyn(k, :, layer) = vpa(subs(stresses_num_yy(:, layer,:), z, i));
%               sigxyn(k, :, layer) = vpa(subs(stresses_num_xy(:, layer,:), z, i));
%               sigxzn(k, :, layer) = vpa(subs(stresses_num_xz(:, layer,:), z, i));
%               sigyzn(k, :, layer) = vpa(subs(stresses_num_yz(:, layer,:), z, i));
%               % sigxxn(k, :, layer) = vpa(subs(stresses_3d_xx(:, :,2), z, i));
%               % sigyyn(k, :, layer) = vpa(subs(stresses_3d_yy(:, :, layer, elem_idx)', z, i));
%               % sigxyn(k, :, layer) = vpa(subs(stresses_3d_xy(:, :, layer, elem_idx)', z, i));
% 
%         depth(k, layer) = i;
%         k = k + 1;
%     end
%     figure (1)
%     plot(sigxxn(:,9,layer), depth(:, layer)/h, "LineStyle","-.","Color","r");
%         xlabel("\sigma_x");
%         ylabel("z*");
%         hold on
% end
% %line 1
% x = [sigxxn(1, 9, 1); sigxxn(end, 9, 2)];
% y = [depth(1, 1)/h; depth(end, 2)/h];
% plot(x, y, 'k');
% %line 1
% x = [sigxxn(1, 9, 2); sigxxn(end, 9, 3)];
% y = [depth(1, 2)/h; depth(end, 3)/h];
% plot(x, y, 'k');
% hold off
