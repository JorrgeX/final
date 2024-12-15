load('mnist.mat'); 

indices_3 = find(labels_train == '3');
indices_8 = find(labels_train == '8');
indices_9 = find(labels_train == '9');
indices = [indices_3; indices_8; indices_9];
labels = labels_train(indices);
X = double(reshape(imgs_train(:, :, indices), [], length(indices))');

m_3 = mean(X(labels == '3', :));
m_8 = mean(X(labels == '8', :));
m_9 = mean(X(labels == '9', :));
m = mean(X);

S_w = zeros(size(X, 2));
X_3 = X(labels == '3', :);
X_8 = X(labels == '8', :);
X_9 = X(labels == '9', :);
S_w = S_w + (X_3 - m_3)' * (X_3 - m_3);
S_w = S_w + (X_8 - m_8)' * (X_8 - m_8);
S_w = S_w + (X_9 - m_9)' * (X_9 - m_9);

n_3 = size(X_3, 1);
n_8 = size(X_8, 1);
n_9 = size(X_9, 1);
S_b = zeros(size(X, 2));
S_b = S_b + n_3 * (m_3 - m)' * (m_3 - m);
S_b = S_b + n_8 * (m_8 - m)' * (m_8 - m);
S_b = S_b + n_9 * (m_9 - m)' * (m_9 - m);

[W, D] = eig(S_b, S_w);
[~, idx] = sort(diag(D), 'descend');
W = W(:, idx(1:2));
Y = X * W;

figure;
scatter(Y(labels == '3', 1), Y(labels == '3', 2), 20, 'r', 'filled');
hold on;
scatter(Y(labels == '8', 1), Y(labels == '8', 2), 20, 'g', 'filled');
scatter(Y(labels == '9', 1), Y(labels == '9', 2), 20, 'b', 'filled');
xlabel('LDA1');
ylabel('LDA2');
legend('3', '8', '9');
title('LDA Projection of Digits 3, 8, and 9');
hold off;

X_centered = X - mean(X);
C = cov(X_centered);
[V_pca, D_pca] = eig(C);
[~, idx_pca] = sort(diag(D_pca), 'descend');
W_pca = V_pca(:, idx_pca(1:2));
Y_pca = X_centered * W_pca;

figure;
scatter(Y_pca(labels == '3', 1), Y_pca(labels == '3', 2), 20, 'r', 'filled');
hold on;
scatter(Y_pca(labels == '8', 1), Y_pca(labels == '8', 2), 20, 'g', 'filled');
scatter(Y_pca(labels == '9', 1), Y_pca(labels == '9', 2), 20, 'b', 'filled');
xlabel('PCA1');
ylabel('PCA2');
legend('3', '8', '9');
title('PCA Projection of Digits 3, 8, and 9');
hold off;