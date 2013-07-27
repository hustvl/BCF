
clear; close all;

load data\mpeg71\1\01_01_feat.mat

i = 20;
x = sx(i,:);
y = sy(i,:);

figure('color', 'w');
subplot(1,2,1), plot( x, y, '.r' );
axis equal

f = feat(:, i);
x = f(1:end/2);
y = f(end/2+1:end);

subplot(1,2,2), plot( x, y, '.r' );
axis equal

% load data\animal_codebook.mat
% 
% for i = 1:100 %size(dict, 2)
%     f = dict(:, i);
%     x = f(1:end/2);
%     y = f(end/2+1:end);
% 
%     figure('color', 'w', 'Position', [1 1 100 100]), plot( x, y, '.-r' );
%     axis equal
%     axis off
%     
%     saveas(gcf, ['fig/codebook/', num2str(i), '.jpg']);
%     close all;
% end