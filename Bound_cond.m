function R=Bound_cond(m,choicebc)

%% Boundary conditions
% Nodes on the boundaries as x0 means @ x=0
for i=1:(2*m+1)
    if i==1
        x0(i)=1;
    else if rem(i,2)==0
        x0(i)=x0(i-1)+(2*m+1);
        else
    x0(i)=x0(i-1)+m+1;
 
end
end
end
x0=x0';
y0=[1:(2*m+1)]';
for i=1:size(x0,1)
    if rem(i,2)==0
    xa(i)=x0(i)+m;
    else
    xa(i)=x0(i)+(2*m);
    end
end
xa=xa';
yb=y0+m*(3*m+2);

rx=[x0;xa];  %Nodes along the x=constant edge
ry=[y0;yb];  %Nodes along the y-constant edge

% % % for simply supported case x0 and xa wil have w and v zero while y0 and yb
% % % will have w and u zero
if choicebc==1
    for q=1:size(rx,1);
        rx1(q)=5*rx(q)-3;
        rx2(q)=5*rx(q)-2;
        rx3(q)=5*rx(q);
    end
    for q=1:size(ry,1);
        ry1(q)=5*ry(q)-4;
        ry2(q)=5*ry(q)-2;
        ry3(q)=5*ry(q)-1;
    end
    R=[rx1 rx2 rx3 ry1 ry2 ry3]';

else if choicebc==2;

        % % %for all edges clamped
        rxy=[x0;xa;y0;yb];
        for q=1:size(rxy,1);
             rxy1(q)=5*rxy(q)-4;
            rxy2(q)=5*rxy(q)-3;
            rxy3(q)=5*rxy(q)-2;
            rxy4(q)=5*rxy(q)-1;
            rxy5(q)=5*rxy(q);
             
        end    
        R=[rxy1 rxy2 rxy3 rxy4 rxy5]';

    else if choicebc==3;

% % for two edges ss two edges x0 and xa and other two free
rxy=[x0;xa];
for q=1:size(rxy,1);
    rxy1(q)=5*rxy(q)-3;
    rxy2(q)=5*rxy(q)-2;
    rxy3(q)=5*rxy(q);
end   
R=[rxy1 rxy2 rxy3]';

        else if choicebc==4;
% % for two edges ss two edges y0 and yb and  other two free
rxy=[y0;yb];
for q=1:size(rxy,1);
     rxy1(q)=5*rxy(q)-3;
    rxy2(q)=5*rxy(q)-2;
    rxy3(q)=5*rxy(q);
end   
R=[rxy1 rxy2 rxy3]';
            end
        end
    end
end