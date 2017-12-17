clear all
close all
clc
%% PCA
N = 20;
mu = [0,0]; % mean
mu2 = [1,1];
sigma = [1,-1.5;1.5,3]; % co-variance
r = mvnrnd(mu,sigma,N); % the data
r2 = mvnrnd(mu2,sigma,N); % the data
c = [ ones(N,1);ones(N,1).*2]; 
X=[r;r2];
c1 = X(find(c==1),:);
c2 = X(find(c==2),:);
figure;

p1 = plot(c1(:,1), c1(:,2), "ro", "markersize",10, "linewidth", 3); hold on;
p2 = plot(c2(:,1), c2(:,2), "go", "markersize",10, "linewidth", 3);

xlim([-5 5]);
ylim([-5 5]);
hold off;
figure;
hold on;
mu = mean(X);
Xm = bsxfun(@minus, X, mu);
C = cov(Xm);
[V,D] = eig(C);
 %sort eigenvectors desc
[D, i] = sort(diag(D), 'descend');
V = V(:,i);
scale = 5;
pc1 = line([mu(1) - scale * V(1,1) mu(1) + scale * V(1,1)], [mu(2) - scale * V(2,1) mu(2) + scale * V(2,1)]);
pc2 = line([mu(1) - scale * V(1,2) mu(1) + scale * V(1,2)], [mu(2) - scale * V(2,2) mu(2) + scale * V(2,2)]);

set(pc1, 'color', [1 0 0], "linestyle", "--");
set(pc2, 'color', [0 1 0], "linestyle", "--");
% project on pc1
z = Xm*V(:,1);
 %and reconstruct it
p = z*V(:,1)';
p = bsxfun(@plus, p, mu);
% delete old plots
%delete(p1);delete(p2);

y1 = p(find(c==1),:);
y2 = p(find(c==2),:);

p1 = plot(y1(:,1),y1(:,2),"ro", "markersize", 10, "linewidth", 3);
p2 = plot(y2(:,1), y2(:,2),"go", "markersize", 10, "linewidth", 3); 
hold off;
figure;
hold on;
% LDA
classes = max(c);
mu_total = mean(X);
mu = [ mean(c1); mean(c2) ];
Sw = (X - mu(c,:))'*(X - mu(c,:));
Sb = (ones(classes,1) * mu_total - mu)' * (ones(classes,1) * mu_total - mu);
[V, D] = eig(Sw\Sb);
% sort eigenvectors desc
[D, i] = sort(diag(D), 'descend');
V = V(:,i);
scale = 5;
pc1 = line([mu_total(1) - scale * V(1,1) mu_total(1) + scale * V(1,1)], [mu_total(2) - scale * V(2,1) mu_total(2) + scale * V(2,1)]);
set(pc1, 'color', [1 0 0], "linestyle", "--");
Xm = bsxfun(@minus, X, mu_total);
z = Xm*V(:,1);
% and reconstruct it
p = z*V(:,1)';
p = bsxfun(@plus, p, mu_total);
% delete old plots
%delete(p1);delete(p2);

y1 = p(find(c==1),:);
y2 = p(find(c==2),:);

p1 = plot(y1(:,1),y1(:,2),"ro", "markersize", 10, "linewidth", 3);
p2 = plot(y2(:,1), y2(:,2),"go", "markersize", 10, "linewidth", 3);
hold off;