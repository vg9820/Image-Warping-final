function IDWImageWarp1(im, psrc, pdst)
%% basic image manipulations
% get image (matrix) size
[h,w,~] = size(im);
% get the number of psrc point
[m,~] = size(psrc);
mu = 2;
im2 = im;
%% Inverse distance-weighted interpolation methods
%% compute the D matrix(2*2*m) of fi
P = psrc;
Q = pdst;
D = zeros(2,2,m);
T = repmat(zeros(2),[m,1]);
% T 2m*2
for i = 1:m
    Pi = repmat(P(i,:),[m,1]);
    Wi = diag(1./power(sum(abs(Pi - P).^2,2).^(1/2),mu));
    % Wi(i,i) is NAN here,replace it into 0
    Wi(i,i) = 0;
    Qi = repmat(Q(i,:),[m,1]);
    Pii = Wi*(P-Pi);
    Qii = Wi*(Q-Qi);
    T((2*i-1):2*i,:) = ((Pii'*Pii)\(Pii'*Qii))';
end   
%% use loop to negate image
%P,Q 2*m matrix
P = P';
Q = Q';
[X1,X2] = meshgrid(1:h,1:w);
p = [X1(:),X2(:)]';

s = zeros(m,h*w);
for i = 1:m
    s(i,:)=1./(sum((p-P(:,i)).^2).^(1/2)).^(mu);
end 
f = zeros(2,h*w);   
for i = 1:m
    wi = s(i,:)./sum(s);
    f = f + wi.*(Q(:,i)+ T((2*i-1):2*i,:)*(p-P(:,i)));
end
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
surface(X2,Y2,Z,'CData', im2,'FaceColor', 'texturemap','Edgecolor','none');




