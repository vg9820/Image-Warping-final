function RBFImageWarp1(im, psrc, pdst)
%% basic image manipulations
% get image (matrix) size
[h, w, ~] = size(im);
% get the number of psrc point
[m, ~] = size(psrc);
% P,Q 2*m
P = psrc';
Q = pdst';
mu =1;
im2 = im;
%% Radial basis functions methods
alpha = zeros(2,m);
f = Q - P;
for i =1:m
    r1(1:m) = sum(abs(P(:,i)-P(:,1:m)).^2).^(1/2);
    r1(i) = 99999;
    r(i) = min(r1);
end
M = zeros(m,m);
for i =1:m
    M(:,i) = (sum(abs(P(:,i)-P).^2)+r.^2).^(mu/2);
end
alpha = (M'\f')'; 
%% use loop to negate image
[X1,X2] = meshgrid(1:h,1:w);
p = [X1(:),X2(:)]';
d = sparse(m,h*w);
for i =1:m
    d(i,:) = sum(abs(P(:,i)-p).^2);
end
f = alpha*(d + r'.^2).^(mu/2) + p;
f = round(f);
f(1,:) = max(f(1,:),1);
f(1,:) = min(f(1,:),h);
f(2,:) = max(f(2,:),1);
f(2,:) = min(f(2,:),w);

X2 = reshape(f(1,:),[h,w]);
Y2 = reshape(f(2,:),[h,w]);
Z = zeros(h,w);
%% get image
ax2 = subplot(122);
cla(ax2);
surface(X2,Y2,Z,'CData', im2,'FaceColor', 'texturemap', 'Edgecolor','none');

