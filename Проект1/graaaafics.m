clc;clear;

M = dlmread('fine.txt');
% t, c 	
% V, м/с 
% tetta, град
% X, м          Y, м        Z, м
% omegaZ, рад/c
% m, kg
% Xc_1, м       Yc_1, м     Zc_1, м
% Xc_2, м       Yc_2, м     Zc_2, м
% n_Ya, -
% fi, град      hi, град    r_viz м
% alfa, град        betta, град
% delta_v, град     delta_n, град

t(:) = M(:,1); 	
V(:) = M(:,2);
tetta(:) = M(:,3);
X(:) = M(:,4);
Y(:) = M(:,5);
Z(:) = M(:,6);
omegaZ(:) = M(:,7);
m(:) = M(:,8);
Xc_1(:) = M(:,9);
Yc_1(:) = M(:,10);
Zc_1(:) = M(:,11);
Xc_2(:) = M(:,12);
Yc_2(:) = M(:,13);
Zc_2(:) = M(:,14);
n_Ya(:) = M(:,15);
fi(:) = M(:,16);
hi(:) = M(:,17);
r_viz(:) = M(:,18);
alfa(:) = M(:,19);
betta(:) = M(:,20);
delta_v(:) = M(:,21);
delta_n(:) = M(:,22);
dr_dt1(:) = M(:,23);
dr_dt2(:) = M(:,24);

% plot(X,Y,Xc_1,Yc_1);
% grid;

% plot(t,delta_v);
% plot(t,delta_n);
% plot(t,dr_dt1-dr_dt2)
% plot(X,Y,Xc_1,Yc_1,Xc_2,Yc_2);
% plot(t,r_viz)

subplot(3,1,1);
plot(t,dr_dt1); grid;
subplot(3,1,2);
plot(t,dr_dt2); grid;
subplot(3,1,3);
plot(t,dr_dt1-dr_dt2); grid;

% plot(t,n_Ya)

% plot(X(160:length(X)),Y(160:length(Y)),Xc_1(160:length(Xc_1)),Yc_1(160:length(Yc_1)));
        % plot(t(1:length(t)-1),diff(fi)./diff(t));

        % plot(t(1:length(t)-1),diff(r_viz)./diff(t));

        %  plot(diff(Xc_1) - V(2:length(t)).*cosd(tetta(2:length(t))), fi(2:length(t)))

% plot(V(161:length(t)).*cosd(tetta(161:length(t))) - diff(Xc_1(160:length(Xc_1))), fi(161:length(t)))
% plot(diff(diff(Xc_1(160:length(Xc_1))) - V(161:length(t)).*cosd(tetta(161:length(t)))), diff(fi(161:length(t))))

% plot(fi(160:length(t)), n_Ya(160:length(t)))

% plot(diff(fi(160:length(t))), diff(n_Ya(160:length(t))))


% plot(t(33:length(t)),diff(r_viz(32:length(t))))
% disp(mean(diff(r_viz(32:length(t)-1))))
