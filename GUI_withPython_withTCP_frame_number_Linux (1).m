close all; clear all; clc;
% Open terminal and run python script


% if python error run this in terminal: pkill -f DataSave2_512.py

% wdir = pwd;
% pythonScript = fullfile(wdir, 'DataSave2_512.py'); %NAME OF FILE TO RUN 
% pythonCmd = sprintf('gnome-terminal -- bash -c "python3 ''%s''; exec bash" &', pythonScript); %WHICH LANGUAGE
% system(pythonCmd);
% disp('Python UDP listener started in background.');

set(0, 'DefaultUIControlFontSize', 12);

%%%%%%%%  init GUI  %%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
set(gcf, 'Position',  [10, 100, 1500, 600]);
c_fre = uicontrol('Style','edit','Position',[1200 555 100 30],'String',2400);
t1 = uicontrol('Style','text','Position',[1030 550 150 40],'String','NCO center Frequency (MHz)');

ROI1 = uicontrol('Style','edit','Position',[1200 485 100 30],'String',680);
ROI2 = uicontrol('Style','edit','Position',[1310 485 100 30],'String',730);

t2 = uicontrol('Style','text','Position',[1010 482 180 40],'String','DF region of interest frequency index');

t3 = uicontrol('Style','text','Position',[180 510 180 30],'String','AZ direction');
t4 = uicontrol('Style','text','Position',[620 510 180 30],'String','EL direction');


CFAR_number = uicontrol('Style','edit','Position',[1430 20 30 20],'String',30);
t5 = uicontrol('Style','text','Position',[1310 20 120 20], 'String','CFAR min limit');

f_number = uicontrol('Style','edit','Position',[1230 20 50 20],'String',512);
t5 = uicontrol('Style','text','Position',[1110 20 120 20], 'String','Frame number');

num_track = uicontrol('Style','edit','Position',[510 20 40 20],'String',7);
t6 = uicontrol('Style','text','Position',[390 20 120 20], 'String','Num of Track');

trail_length = str2num(num_track.String);
trail = NaN(trail_length, 2, 2); %trail length, XY, AZ/EL
color_int = (trail_length:-1:1)/trail_length;
C1 = color_int(:);
C_norm = (C1 - min(C1)) / (max(C1) - min(C1));
RGB = [C_norm, zeros(length(C_norm),1), zeros(length(C_norm),1)];


disp('here')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


frequency_lib = 2456; %MHz (use 2456 MHz or 2360 MHz calibration library)
if(frequency_lib==2456)
    libpath = './Calibration_Library_2p45GHz';
elseif(frequency_lib==2360)
    libpath = './Calibration_Library_2p36GHz';
end

% ================= LOAD CALIBRATION LIBRARY =================
lib_cache = fullfile(libpath, 'cached_library_data.mat');
if isfile(lib_cache)
    load(lib_cache, 'Lib_Mag', 'Lib_Phase', 'Lib_Complex');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
algorithm_select='CAPON';

plot_flag=0;      % Enable plotting
print_detail=0;   % Print debug/info
save_flag=0;      % Save results


num_frame = 512;
N_frame=num_frame;     % Number of frames
center_fre = str2num(c_fre.String);

serverIP = '192.168.56.10';
serverPort = 9000;
% Try to establish a single persistent TCP connection
while true
    try
        t = tcpclient(serverIP, serverPort);
        disp('Connected to TCP server.');
        break; % Exit loop once connected
    catch ME
        warning('Connection failed: %s\nRetrying in 2 seconds...', ME.message);
        pause(2);
    end
end


% ================= MAIN LOOP =================
for ang_ind=1:10000   % Placeholder loop (expandable if scanning angles)
tic   

    fprintf('--- Receiving frame #%d ---\n', ang_ind);

    % Preallocate buffer for all data
    % dataAll = int16.empty(totalSize, 0);


%%%%%%%%%%%%%%%%%%%%   file load   %%%%%%%%%%%%%%%%%%%%%%%

% ang_ind=18;
fileID = fopen(['Falcon_Raw/' num2str(ang_ind,'%04d') '.bin'], 'r', 'ieee-le');


time_out=0;
while (fileID<0)&(time_out<100)
    pause(0.05);
    fileID = fopen(['Falcon_Raw/' num2str(ang_ind,'%04d') '.bin'], 'r', 'ieee-le');
    time_out=time_out+1;
end

if (time_out<100)
    pause(0.05);
    C = fread(fileID, Inf, 'int16');

    time_out2=0;
    while (length(C)~=8388608/512*num_frame )&(time_out2<20)
        pause(0.05)
         fileID = fopen(['Falcon_Raw/' num2str(ang_ind,'%04d') '.bin'], 'r', 'ieee-le');
        C = fread(fileID, Inf, 'int16');
        time_out2=time_out2+1;
    end
    if (time_out2>=20)
        fprintf('Time out, file size not correct \n');
    end
    fclose(fileID);
else
    fprintf('Time out, no file received \n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_all=[];


C0 = reshape(C,[2,length(C)/2]).';
L_C0=length(C0);

C1= C0 (1:2:L_C0/4,:).'; C_all(:,1)=C1(:);
C2= C0 (2:2:L_C0/4,:).'; C_all(:,2)=C2(:);
C3= C0 (L_C0/4+1:2:L_C0/4*2,:).'; C_all(:,3)=C3(:);
C4= C0 (L_C0/4+2:2:L_C0/4*2,:).'; C_all(:,4)=C4(:);
C5= C0 (L_C0/2+1:2:L_C0/4*3,:).'; C_all(:,5)=C5(:);
C6= C0 (L_C0/2+2:2:L_C0/4*3,:).'; C_all(:,6)=C6(:);
C7= C0 (L_C0/4*3+1:2:L_C0,:).'; C_all(:,7)=C7(:);
C8= C0 (L_C0/4*3+2:2:L_C0,:).'; C_all(:,8)=C8(:); 

C1_cmplex=C_all(1:2:end,:)+1i*C_all(2:2:end,:);

C1_cmplex_ooo=C1_cmplex;

C1_cmplex=C1_cmplex(length(C1_cmplex)-num_frame*1024+1:length(C1_cmplex),:);

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fre_sample=2.94912e9;        % Original sample rate
    fre_sample=2.94912e9/12;     % Decimated sample rate
    fre=[0:1023]/1024*fre_sample; % Frequency vector

    roi_1=str2num(ROI1.String);
    roi_2=str2num(ROI2.String);
    roi_freq=[roi_1:roi_2];

if trail_length~=str2num(num_track.String)
    trail_length = str2num(num_track.String);
    trail = NaN(trail_length, 2, 2); %trail length, XY, AZ/EL
    color_int = (trail_length:-1:1)/trail_length;
    C = color_int(:);
    C_norm = (C - min(C)) / (max(C) - min(C));
    RGB = [C_norm, zeros(length(C_norm),1), zeros(length(C_norm),1)];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ================= FFT PROCESSING =================
    FFT_full=zeros(N_frame,1024,6);      % FFT results [frame, freq, RX]
    threshold_FFT=zeros(N_frame,1024,6); % Thresholding results
   
window_hann=hann(1024);



    for frame_ind=1:N_frame
        % Extract FFT per channel (mapping order based on antenna channels)
        % freA=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),1));
        freB=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),2)); % channel 5
        freC=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),3)); % channel 4
        freD=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),4)); % channel 6
        freE=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),5)); % channel 1
        freF=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),6)); % channel 2
        freG=fft(window_hann.*C1_cmplex([1:1024]+1024*(frame_ind-1),7)); % channel 3
        % freH=fft(C1_cmplex([1:1024]+1024*(frame_ind-1),8));
        
        % Reorder into 6 RX channels
        FFT_full(frame_ind,:,1)=fftshift(freE);
        FFT_full(frame_ind,:,2)=fftshift(freF);
        FFT_full(frame_ind,:,3)=fftshift(freG);
        FFT_full(frame_ind,:,4)=fftshift(freC);
        FFT_full(frame_ind,:,5)=fftshift(freB);
        FFT_full(frame_ind,:,6)=fftshift(freD);      
    end


    % ================= PLOT SPECTROGRAM =================
    
    
    Nrx=6;   % Number of RX channels used
    time_axis=[0:1023]*1/fre_sample*1e6;  % Time axis (us)
    % fre_axis=fre/1e9+2.277;               % Frequency axis (GHz), with LO offset
    fre_axis=fre/1e9-fre_sample/1e9/2;               % Frequency axis (GHz), with LO offset
if(1)    
    for rx=3:3
        test_STFT=abs(FFT_full(:,:,rx));  % Magnitude spectrogram for RX
        
        
            figure(1);
            % subplot(2, 3, rx);
            subplot('Position',[0.68 0.2 0.3 0.5]);

            imagesc(time_axis,fre_axis,mag2db(test_STFT.')); % Plot spectrogram
            set(gca, 'YDir', 'normal');  % Flip Y axis
            hold on;
            plot(time_axis,ones(1,1024)*fre_axis(roi_1),'r:','linewidth',3);
            plot(time_axis,ones(1,1024)*fre_axis(roi_2),'r:','linewidth',3);
            hold off;

            grid on;
            
            xlabel('Time (us)');
            ylabel('Frequency (GHz)');
            cb = colorbar;               
            ylabel(cb, 'Magnitude (dB)');    
            caxis([20 90]);              % Dynamic range for plotting
            title(sprintf('Spectrogram of RX %d', rx));
    end
end
     
 
    offset=1.25;

    mag=abs(FFT_full);
    noise_all0=median(mag,2);
    noise_all=repmat(noise_all0,[1 1024 1]);

    threshold_cfar = offset * noise_all;
    threshold_FFT = mag > threshold_cfar;


    % ================= COMBINE DETECTIONS ACROSS RX =================
    final_threshold=ones(N_frame,1024);
    for rx=1:Nrx
        test_threshold=threshold_FFT(:,:,rx);
        
        if(plot_flag)
            figure(97);
            subplot(2, 3, rx);
            imagesc(time_axis,fre_axis,(test_threshold.')); % Binary detections
            set(gca, 'YDir', 'normal');
            grid on;
            xlabel('Time (us)');
            ylabel('Frequency (GHz)');
            colormap(flipud(gray))
            cb = colorbar;              
            title(sprintf('Binary threshold of RX %d', rx));
        end
        
        final_threshold=final_threshold.*test_threshold; % Common detections
    end
    
    % Final binary threshold (intersection across RXs)
    if(plot_flag)
        figure(200)
        imagesc(time_axis,fre_axis,(final_threshold.'));
        set(gca, 'YDir', 'normal');
        grid on;
        xlabel('Time (us)');
        ylabel('Frequency (GHz)');
        colormap(flipud(gray))
        cb = colorbar;              
        title(sprintf('Final threshold'));
    end


    % ================= ANGLE SCAN GRID =================
    if(frequency_lib == 2456)
        AZ_start = 180; AZ_end = -180; AZ_step = -3;
        EL_start = 66; EL_end = 0; EL_step = -3;
    elseif(frequency_lib ==2360)
        AZ_start = 180; AZ_end = -180; AZ_step = -3;
        EL_start = 66; EL_end = -6; EL_step = -3;
    end

    % Build azimuth and elevation scan tables
    AZ_data = AZ_start:AZ_step:AZ_end;
    AZ_steps = length(AZ_data);
    EL_data = EL_start:EL_step:EL_end;
    EL_steps = length(EL_data);
    
    AZ_table=AZ_data;
    EL_table=EL_data;
    length_el=length(EL_table);
    length_az=length(AZ_table);
    
    

    
    % ================= INITIALIZE RESULT ARRAYS =================
    AZ_result=[];
    AZ_result_itp=[];
    EL_result=[];
    EL_result_itp=[];
    summary_result=[];

    % ================= MAIN PEAK PROCESSING LOOP =================
    for peak_idx = 1:length(roi_freq)
        peak = roi_freq(peak_idx);
        fre_peak=fre(peak)/1e9+2.277; % Convert bin index → GHz

        % Find frames where CFAR detected the peak
        frames = find(final_threshold(:,peak) ==1);
        
        if(print_detail)
            fprintf('Peak %d is found in frame(s): %s\n', peak, mat2str(frames'));
        end
        
        if (length(frames)>10)
            % Extract complex FFT data for selected frames
            test_data=squeeze(FFT_full(frames,peak,:));
            X=test_data.'; 
            R = (X * X') / length(frames); % Covariance matrix
            
            % ================= DOA ALGORITHM =================
            if strcmp(algorithm_select,'MUSIC')
                % --------- MUSIC Algorithm ---------
                [Evecs, Evals] = eig(R); % Eigen decomposition
                [Evals_sorted, idx] = sort(diag(Evals), 'descend');
                Evecs_sorted = Evecs(:, idx);
                d=2; % number of signals assumed
                En = Evecs_sorted(:, d+1:end); % Noise subspace
                
                % MUSIC spectrum evaluation
                Pmusic = zeros(length_el, length_az);
                for el_idx = 1:length_el
                    for az_idx = 1:length_az
                        steering = Lib_Complex(:, az_idx, el_idx);  
                        Pmusic(el_idx, az_idx) = 1 / abs(steering' * (En * En') * steering);
                    end
                end
                
                Pmusic_dB = mag2db(abs(Pmusic));
                abs_spectrum=abs(Pmusic);
    
                % Plot MUSIC spectrum if only 1 ROI frequency
                if(length(roi_freq)==1)
                    figure;
                    imagesc(AZ_table, EL_table, Pmusic_dB);
                    set(gca, 'YDir', 'normal');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    grid on; title('2D MUSIC Spectrum');
                    
                    [AZ_grid, EL_grid] = meshgrid(AZ_table, EL_table);
                    figure; surf(AZ_grid, EL_grid, Pmusic_dB, 'EdgeColor', 'none');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Spatial Spectrum (dB)');
                    title('2D MUSIC Spectrum (3D View)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    view(45, 45); shading interp;
                end
                
            elseif strcmp(algorithm_select,'CAPON')
                % --------- CAPON Algorithm ---------
                try
                    R_inv = pinv(R);  % Pseudo-inverse covariance matrix
                catch
                    fprintf('ERROR: Singular matrix\n');
                    ADSINR = -3;  
                    return;
                end
               
                PAD = zeros(length_el, length_az);
                for el_idx = 1:length_el
                    for az_idx = 1:length_az
                        steering = Lib_Complex(:, az_idx, el_idx);
                        PAD(el_idx, az_idx) = (steering' * R_inv) * steering; 
                    end
                end
                PCapon=1./PAD;
                
                abs_spectrum=abs(PCapon);
                PCapon_dB = mag2db(abs(PCapon));
                
                if(length(roi_freq)==1)
                    figure;
                    imagesc(AZ_table, EL_table, PCapon_dB);
                    set(gca, 'YDir', 'normal');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    grid on; title('2D CAPON Spectrum');
                    
                    [AZ_grid, EL_grid] = meshgrid(AZ_table, EL_table);
                    figure; surf(AZ_grid, EL_grid, PCapon_dB, 'EdgeColor', 'none');
                    xlabel('Azimuth (°)'); ylabel('Elevation (°)'); zlabel('Spatial Spectrum (dB)');
                    title('2D CAPON Spectrum (3D View)');
                    cb = colorbar; ylabel(cb, 'Spatial Spectrum (dB)');
                    view(45, 45); shading interp;
                end
                
            elseif strcmp(algorithm_select,'Correlation')
                % --------- Correlation-based DOA ---------
                rxSignal2D_FFT_rx1=test_data(:,1);
                rxSignal2D_FFT_test=test_data./rxSignal2D_FFT_rx1;
                rxSignal2D_FFT_test_mean=mean(rxSignal2D_FFT_test,1).';
                weighting=[1 1 1 1 1 1]';
                
                for az=1:AZ_steps
                    for el=1:EL_steps
                        r=corrcoef(rxSignal2D_FFT_test_mean.*weighting, ...
                                   squeeze(Lib_Complex(:,az,el)).*weighting);
                        coe(az,el)=r(1,2);
                    end
                end
                
                abs_spectrum=reshape(abs(coe),[length_az length_el]);
                
                if(length(roi_freq)==1)
                    figure();
                    imagesc(EL_table, AZ_table, mag2db(abs_spectrum));
                    set(gca, 'YDir', 'normal');
                    xlabel('Elevation (deg)'); ylabel('Azimuth (deg)');
                    colorbar; colormap jet;
                    title('Color Field Map');
                    grid on;
                end
                abs_spectrum=abs_spectrum.';
            end
           
            % ================= PEAK SEARCH & INTERPOLATION =================
            [max_val, max_idx] = max(abs_spectrum(:));
            [peak_el, peak_az] = ind2sub(size(abs_spectrum), max_idx);
            
            AZ_peak = AZ_table(peak_az);
            EL_peak = EL_table(peak_el);
            
            cf_reshape = abs_spectrum;       
            AZ_step = AZ_table(2) - AZ_table(1);
            EL_step = EL_table(2) - EL_table(1);
            
            % Interpolation in Azimuth
            delta_az = 0;
            if peak_az > 1 && peak_az < length_az
                a = cf_reshape(peak_el, peak_az - 1);
                b = cf_reshape(peak_el, peak_az);
                c = cf_reshape(peak_el, peak_az + 1);
                if a ~= c
                    delta_az = 0.5 * (a - c) / (a - 2*b + c);
                end
            end
            
            % Interpolation in Elevation
            delta_el = 0;
            if peak_el > 1 && peak_el < length_el
                a1 = cf_reshape(peak_el - 1, peak_az);
                b1 = cf_reshape(peak_el, peak_az);
                c1 = cf_reshape(peak_el + 1, peak_az);
                if a1 ~= c1
                    delta_el = 0.5 * (a1 - c1) / (a1 - 2*b1 + c1);
                end
            end
            
            % Final interpolated peaks
            AZ_peak_itp = AZ_table(peak_az) + delta_az * AZ_step;
            EL_peak_itp = EL_table(peak_el) + delta_el * EL_step;
                        
            % Store results
            AZ_result=[AZ_result;AZ_peak];
            AZ_result_itp=[AZ_result_itp;AZ_peak_itp];
            EL_result=[EL_result;EL_peak];
            EL_result_itp=[EL_result_itp;EL_peak_itp];
        
            summary_result=[summary_result;peak fre_peak length(frames) AZ_peak EL_peak AZ_peak_itp EL_peak_itp];

        end
    end


    if ~isempty(AZ_result)
        final_AZ=median(AZ_result);
        final_AZ_ITP=median(AZ_result_itp);
    
        final_EL=median(EL_result);
        final_EL_ITP=median(EL_result_itp);
        if(print_detail)
            fprintf('Final Results:\n');
            fprintf('  AZ      = %.3f\n', final_AZ);
            fprintf('  AZ ITP  = %.3f\n', final_AZ_ITP);
            fprintf('  EL      = %.3f\n', final_EL);
            fprintf('  EL ITP  = %.3f\n', final_EL_ITP); 
        end
        if(save_flag==1)
            if(frequency_lib==2456)
                save_folder='Result_2p45_lib';
                save_folder = [save_folder '_' algorithm_select];
                if ~exist(save_folder, 'dir')
                    mkdir(save_folder);
                end
                save([save_folder '/' num2str(id_name) '_2p45_result.mat'], 'summary_result','final_AZ','final_AZ_ITP','final_EL','final_EL_ITP')
                id_name=id_name+1;
            elseif(frequency_lib==2360)
                save_folder='Result_2p36_lib';
                save_folder = [save_folder '_' algorithm_select];
                if ~exist(save_folder, 'dir')
                    mkdir(save_folder);
                end
                save([save_folder '/' num2str(id_name) '_2p36_result.mat'], 'summary_result','final_AZ','final_AZ_ITP','final_EL','final_EL_ITP')
                id_name=id_name+1;
            end
        end
    end

    if ~isempty(AZ_result)
    distance=1;
    X_plot=distance*cosd(final_AZ_ITP+90);Y_plot=distance*sind(final_AZ_ITP+90);
    X_plot2=distance*cosd(final_EL_ITP);Y_plot2=distance*sind(final_EL_ITP);

        
        if(length(summary_result(:,3))>str2num(CFAR_number.String))    
            
            for i = trail_length:-1:2
                trail(i, :, 1) = trail(i-1, :, 1);
                trail(i, :, 2) = trail(i-1, :, 2);
            end
            trail(1, :, 1) = [X_plot, Y_plot];
            trail(1, :, 2) = [X_plot2, Y_plot2];
            
            figure(1)
            subplot('Position',[0.05 0.2 0.25 0.63]);

            % hold on;
            plot(X_plot,Y_plot,'xk',[0,X_plot],[0,Y_plot],'r:','MarkerSize',12,'Linewidth',2);
            hold on
            % scatter(squeeze(trail(:, 1, 1)),squeeze(trail(:, 2, 1)), 36, RGB, 'filled', 'MarkerFaceAlpha', 'flat', 'AlphaData', C_norm);
            for i = 1:trail_length - 1
                if ~(isnan(trail(i:i+1, :, 1)))
                    plot([trail(i, 1, 1) trail(i+1, 1, 1)],[trail(i, 2, 1) trail(i+1, 2, 1)],'Color',[1 0 0 C_norm(i)],'LineWidth',2);
                end
            end
            hold off
            grid on;
            xlim([-1.5 1.5]);
            ylim([-1.5 1.5]);
            xlabel('X');ylabel('Y');
            text(-0.5, 1.35, ['Azimuth Angle = ' num2str(-round(90 - (final_AZ_ITP + 90), 2))], 'FontSize', 15, 'Color', 'r')

            
            subplot('Position',[0.35 0.2 0.25 0.63]);

            % hold on;
            plot(X_plot2,Y_plot2,'xk',[0,X_plot2],[0,Y_plot2],'r:','MarkerSize',12,'Linewidth',2);
            hold on
            % scatter(squeeze(trail(:, 1, 2)),squeeze(trail(:, 2, 2)), 36, RGB, 'filled', 'MarkerFaceAlpha', 'flat', 'AlphaData', C_norm);
            for i = 1:trail_length - 1
                if ~(isnan(trail(i:i+1, :, 2)))
                    plot([trail(i, 1, 2) trail(i+1, 1, 2)],[trail(i, 2, 2) trail(i+1, 2, 2)],'Color',[1 0 0 C_norm(i)],'LineWidth',2);
                end
            end            
            hold off
            grid on;
            xlim([-0.5 1.5]);
            ylim([-0.5 1.5]);
            xlabel('X');ylabel('Z');
            text(0.3, 1.4, ['Elevation Angle = ' num2str(round(final_EL_ITP, 2))], 'FontSize', 15, 'Color', 'r')

            % fontsize(12, "points");
        else
            fprintf('No drone signal \n');
         end    
    else
    
        fprintf('No detection:\n');
    end


% 
% 
% 
% 
% %%%%%%%%%%%%%%%%   add python capture next file  %%%%%%%%%%%%%%%%
% Freq = [2400 2200 2000 1800 1600 1300 900];
% 
if center_fre ~= str2num(c_fre.String)

 cmd = sprintf('SetSingleFrequency %d\n', str2num(c_fre.String));
 write(t, uint8(cmd), "uint8");  % Send as raw bytes
 center_fre = str2num(c_fre.String);
    fprintf('NCO updated \n');
end
 %pause(0.03);
% 
% %cmd = sprintf('SetRangeFrequency %d %d\n',1,4);
% %write(t, uint8(cmd), "uint8");  % Send as raw bytes
% 
% 

read_frame_number = str2num(f_number.String);
 if  num_frame ~= read_frame_number
    if ismember(read_frame_number, [128 256 512 1024])  
         cmd = sprintf('SetDwellSize %d\n', read_frame_number);
         write(t, uint8(cmd), "uint8");  % Send as raw bytes
         num_frame = read_frame_number;
         N_frame=num_frame;
    end

 end
 % 
 %  if  ang_ind == 100
 % cmd = sprintf('SetDwellSize %d\n',512);
 % write(t, uint8(cmd), "uint8");  % Send as raw bytes
 % num_frame = 512;
 % end
 


 cmd = sprintf('NextFrame\n');
 write(t, uint8(cmd), "uint8");  % Send as raw bytes
% 
% %write(t, uint8('NextFrame\n'), "uint8");  % Send as raw bytes
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc

end
