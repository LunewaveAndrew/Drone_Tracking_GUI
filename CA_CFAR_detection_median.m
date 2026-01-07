% num_training = 70;          % Number of training cells on each side
% num_guard = 5;              % Number of guard cells on each side
% offset = 2.7;                 % Threshold scaling factor Good

offset = 10;                 % Threshold scaling factor Good

% Initialize threshold and detection vectors
threshold_cfar = 0;
cfar_output = zeros(1, length(mag));


% Apply CA-CFAR
noise_level = median(mag);
threshold_cfar = offset * noise_level;
cfar_output = mag > threshold_cfar;
% 
% for i = num_training + num_guard + 1 : length(mag) - (num_training + num_guard)
%     % Define training cells
%     start1 = i - num_guard - num_training;
%     end1   = i - num_guard - 1;
%     start2 = i + num_guard + 1;
%     end2   = i + num_guard + num_training;
% 
%     training_cells = [mag(start1:end1), mag(start2:end2)];
% 
%     % Estimate noise level
%     noise_level = mean(training_cells);
% 
%     % Calculate threshold
%     threshold_cfar(i) = offset * noise_level;
% 
%     % Compare CUT to threshold
%     if mag(i) > threshold_cfar(i)
%         cfar_output(i) = 1;
%     end
% end


