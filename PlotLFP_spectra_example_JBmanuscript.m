% plot example LFP and behavior - 1 hr block
% January 2023



% LOAD JB26 LFPinfo
% LOAD JB26 STATETIMES
%load('JB26_LFPinfo.mat')
%load('JB26_STATETIMES.mat')
%statetimes = STATETIMES;

new_statetimes = statetimes;
 % get rid of duplicate states
     testD   = diff(new_statetimes(:,1));
     killem  = find(testD == 0);
     killem  = killem+1;
     new_statetimes(killem,:) = [];

    % find and eliminate transient wake episodes (e.g. the animal wakes up
    % simply to shift position)
    nsec=60; %set limit at 1 min-eliminate all states that only last 1 minute or less
    sdur = diff(new_statetimes(:,2));
    fcks = find(new_statetimes(1:end-1,1)>2 & sdur<=nsec ); % find instances of waking that are less than nsec seconds long

    tempf = new_statetimes(fcks+1,1)<3; % find one fcks for which the following state is still sleeping
    fcks2 = fcks(tempf);
    new_statetimes(fcks2,:) = [];
    % second pass
    testD   = diff(new_statetimes(:,1));
    killem  = find(testD == 0);
    killem  = killem+1;
    new_statetimes(killem,:) = [];
    
    

  i = 80; % choose block of interest
  total_powers = sum(LFPinfo(i).spectrogram.P);
  max_powers = max(LFPinfo(i).spectrogram.P);
  powerinFreqs = LFPinfo(i).spectrogram.P ./ max_powers;
  
   figure()
    imagesc([0:.016:60],(LFPinfo(i).spectrogram.F),10*log10(abs(LFPinfo(i).spectrogram.P)));
    %imagesc([0:.016:60],(LFPinfo(i).spectrogram.F),10*log10(powerinFreqs));
    axis([0 60 0 15]);
    colorbar

    ylabel ('Frequency')
    xlabel ('Time (minutes)')
    title (i)
    set(gca,'YDir','normal', 'FontSize', 15)
    
    
LFPstart = LFPinfo(i).startTime;
first_state = find(new_statetimes(:,2)>LFPstart, 1);
last_state = find(new_statetimes(:,2)<(LFPstart+3600), 1, 'last');

figure()
hold on
xs = [0:.016:60];
these_statetimes = new_statetimes(first_state:last_state, :);
these_statetimes(:,2) = (these_statetimes(:,2) - LFPstart)./60;
for i = 1:length(these_statetimes(:,1))
    
    if these_statetimes(i,1) == 1 % REM
        FC = 'r';
    elseif these_statetimes(i,1) == 2 % NREM
        FC = 'b';
    elseif these_statetimes(i, 1) == 5 % QW
        FC = 'g';
    elseif these_statetimes(i,1) == 4 % AW
        FC = 'y';
    end
    
    if i~= length(these_statetimes(:,1))
        rectangle('Position', [these_statetimes(i,2), 1, (these_statetimes(i+1,2)-these_statetimes(i,2)), 2], 'FaceColor', FC)
    else
        rectangle('Position', [these_statetimes(i,2), 1, (3600-these_statetimes(i,2)), 2], 'FaceColor', FC)
    end

end
xlim([0,60])

