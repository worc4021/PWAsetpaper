clear 
close all

n = 10;
R = linspace(-2,2,n);
Y = R.^3*.1;

d = (R(end)-R(1))/(n-1);

m = diff(Y)/d;
b = Y(2:end)-m.*R(2:end);

figure(1)
plot(R,R.^3*.1);
hold('on');
fplot(@(x).1*x.^3,[R(1),R(end)],'linewidth',2);
for i = 1:5
    plot([-2,0],m(i)*[-2,0]+b(i))
end
for i = 5:9
    plot([0,2],m(i)*[0,2]+b(i))
end
plot([-2,0],[0,0],[0,2],[0,0])
hold off
ylim([-1,1])
xlabel('$t$')
ylabel('$f_1(t)$')

% Z = exp(-R/4)-1;
% 
% M = zeros(size(Z));
% M(1:5) = diff(Z(1:6))/d;
% B = zeros(size(Z));
% B(1:5) = Z(2:6)-M(1:5).*R(2:6);
% M(6:end) = -1/4*exp(-R(6:end)/4);
% B(6:end) = Z(6:end)-M(6:end).*R(6:end);
% 
% 
% figure(2)
% fplot(@(x)exp(-x/4)-1,[-2,2],'linewidth',2)
% hold on
% for i = 1:5
%     plot([-2,0],M(i)*[-2,0]+B(i))
% end
% for i = 6:10
%     plot([0,2],M(i)*[0,2]+B(i))
% end
% hold off

heav = @(t)(1+sign(t))/2;
be = @(t)1/2*asin(sqrt(abs(t)/2))-t;

R = [-2,-1.8,-1.5,0,.7,1,1.5,1.7,2];

Z = real(1/2*asin(heav(R).*sqrt(abs(R)/2) + heav(-R).*R/2));
M = diff(Z)./diff(R);
B = Z(1:end-1)-M.*R(1:end-1);
% M(1:3) = diff(Z(1:4))./diff(R(1:4));
% B(1:3) = Z(1:3)-M(1:3).*R(1:3);
% M(4:9) = 1./(4*sqrt(-(R(4:end-1)+diff(R(4:end))/2-4).*(R(4:end-1)+diff(R(4:end))/2)));
% B(4:9) = 1/2*asin(sqrt(abs(R(4:end-1)+diff(R(4:end))/2)/2) )-(R(4:end-1)+diff(R(4:end))/2).*M(4:end);

figure(2)
fplot(@(x)real(1/2*asin(heav(x).*sqrt(abs(x)/2) + heav(-x).*x/2)),[-2,2],'linewidth',2)
hold('on');
for i = 1:3
    plot([-2,0],M(i)*[-2,0]+B(i))
end
for i = 5:8
    plot([0,2],M(i)*[0,2]+B(i))
end
plot([-2,0],[0,0],[0,2],[0,0])
hold off
ylim([-1,1])
xlabel('$t$')
ylabel('$f_2(t)$')

M1 = [0,0;kron(m(5:end)',[1/2,1])];
B1 = [0;b(5:end)'];
M2 = [0,0;kron(M(5:end)',[1,-1/2])];
B2 = [0;B(5:end)'];
M3 = [0,0;kron(-m(1:5)',[1/2,1])];
B3 = [0;-b(1:5)'];
M4 = [0,0;kron(-M(1:3)',[1,-1/2])];
B4 = [0;-B(1:3)'];

Q = -blkdiag(ones(size(B1)),ones(size(B2)),ones(size(B3)),ones(size(B4)));

T = [1,0,0,0;
    0,1,0,0;
    1,0,0,0;
    0,0,0,-1;
    0,0,-1,0;
    0,1,0,0;
    0,0,-1,0;
    0,0,0,-1];
L = [1,1;-1,1;1,-1;-1,-1];
l = ones(4,1)*2;

D = zeros(4,4,4);
D(:,:,1) = kron(eye(4),L(1,:))*T;
D(:,:,2) = kron(eye(4),L(2,:))*T;
D(:,:,3) = kron(eye(4),L(3,:))*T;
D(:,:,4) = kron(eye(4),L(4,:))*T;

% A = [.2,.3;-.3,.2];


    Aineq = [zeros(16,2),kron(eye(4),-ones(4,1)),[D(:,:,1);D(:,:,2);D(:,:,3);D(:,:,4)];...
            [M1;M2;M3;M4],zeros(size(Q,1),4),Q;...
            L,diag([1,1,1,1]),zeros(4)];

    bineq = [zeros(16,1);-[B1;B2;B3;B4];l];

    [VAU,t] = LRS(struct('rep','H','Aineq',Aineq,'bineq',bineq));

    VV = bigVReduce(VAU(logical(t),1:2));
    vv = zeros(size(VV,1)*4,2);
    fv = zeros(size(VV));
    for i = 1:size(VV,1)
        x = VV(i,:)';
        t = [max(M1*x+B1),max(M2*x+B2),max(M3*x+B3),max(M4*x+B4)]';
        ValidateF(x,t);
        dx = reshape(T*t,2,4)+repmat(x,1,4);
        vv((4*i-3):4*i,:) = dx';
        fprintf('\\draw[cyan] (%8.4f,%8.4f) -- (%8.4f,%8.4f) -- (%8.4f,%8.4f) -- (%8.4f,%8.4f) -- cycle;\n',...
            dx(1),dx(2),dx(3),dx(4),dx(7),dx(8),dx(5),dx(6));
        r = [1,-1/2]*x;
        dx = [([1/2,1]*x)^3/10;
            real(1/2*asin(heav(r).*sqrt(abs(r)/2) + heav(-r).*r/2))];
        fv(i,:) = VV(i,:)+dx';
    end
    
    k = convhull(VV(:,1),VV(:,2));
    s = [1/2,1;1,-1/2];
    achs = s*[eye(2),-eye(2)]*(sqrt(2)*1.1);
    
    figure(10)
    plot(VV(k,1),VV(k,2),...
        VV(:,1),VV(:,2),'x',...
        vv(:,1),vv(:,2),'o',...
        2*[1,0,-1,0,1],[0,1,0,-1,0]*2,...
        fv(:,1),fv(:,2),'x',...
        achs(1,[1,3]),achs(2,[1,3]),'-.',...
        achs(1,[2,4]),achs(2,[2,4]),'-.')
xlabel('$x_1$');
ylabel('$x_2$');

% fprintf('\n');
% for i = 1:size(vv,1)
%     fprintf('%8.4f/%8.4f,',vv(i,1),vv(i,2));
% end
% fprintf('\n');
Vn = [];
for i = 1:4
    IDX = false(4,1);
    IDX(i) = true;
    Aineq = [zeros(4,2),-ones(4,1),D(:,:,i);...
            [M1;M2;M3;M4],zeros(size(Q,1),1),Q;...
            kron(L,ones(4,1)),zeros(16,1),[D(:,:,1);D(:,:,2);D(:,:,3);D(:,:,4)]];
	bineq = [zeros(4,1);-B1;-B2;-B3;-B4;kron(l,ones(4,1))];
    Aeq = [L(IDX,:),1,zeros(1,4)];
    beq = l(IDX);
    [V,t] = LRS(struct('rep','H','Aineq',Aineq,'bineq',bineq,'Aeq',Aeq,'beq',beq));
    Vn = [Vn;bigVReduce(V(logical(t),1:2))];
end