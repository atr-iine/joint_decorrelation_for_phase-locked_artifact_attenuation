function [outEEG, errorflag] = jointdecorrelation(inEEG, thresh)
% Joint correlation
%   [outEEG, errorflag] = jointdecorrelation(inEEG, thresh)
%
%   Inputs:
%       inEEG: EEG struct of EEGLAB
%       thresh: Threshold for BCG correction
%   Outputs:
%       outEEG: EEG struct of EEGLAB
%       errorflag: If true, some error has occurred in computation; otherwise, no error
%
%   The original Python program: (c) Reinmar J. Kobler, 2024
%   This MATLAB version: (c) Toshikazu Kuroda, 2024
%
%
%   Reference: Kuroda, T., Kobler, R.J., Ogawa, T., Tsutsumi, M.,
%   Kishi, T., & Kawanabe, M. (2024). Test-retest reliability of EEG
%   microstate metrics for evaluating noise reductions in simultaneous
%   EEG-fMRI. Imaging Neuroscience. https://doi.org/10.1162/imag_a_00272 
%

    try
        ecg_evt_ixs = find(ismember({inEEG.event.type},{'fmrib_rpeak'}));
        tmp = struct2cell(inEEG.event);
        latency = cell2mat(squeeze(tmp(2,1,:)));
        ecg_trig = int64(round(latency(ecg_evt_ixs))); % Note: Latency has already been adjusted in our case
    
        offset = int64(-0.2*inEEG.srate);
        duration = int64(0.7*inEEG.srate);    
        I = ecg_trig + (0:(duration-1)) + offset;
        I = I(all((I > 0) & (I < size(inEEG.times,2)),2),:);

   
        rnd_bins = int64(randi([-offset,(latency(end) - (duration + offset))], size(ecg_trig,1), size(ecg_trig,2)));          
        IR = rnd_bins + (0:(duration-1)) + offset; 
        IR = IR(all((IR > 0) & (IR < size(inEEG.times,2)),2),:); 
        IR = IR(1:size(I,1),:);
    
        x_art = art_ref(inEEG,I);
        x_ref = art_ref(inEEG,IR);
    
        x_ref_flat = reshape(x_ref,size(inEEG.data,1),[],1);
        cov_x = (x_ref_flat * x_ref_flat') / size(x_ref_flat,2);
        
        n_epo = size(x_art,3);
        n_epo_avg = 25; 
        n_obs_avg = floor(n_epo/n_epo_avg);
        n_epo = n_obs_avg*n_epo_avg;

        tmp0 = x_art(:,:,1:n_epo);
        tmp1 = reshape(tmp0,[size(x_art,1),size(x_art,2), n_epo_avg, n_obs_avg]); 
        tmp2 = permute(tmp1, [1,2,4,3]);
        tmp3 = mean(tmp2,4);
        tmp4 = permute(tmp3,[1,3,2]);
        x_art = reshape(tmp4,size(tmp4,1),[],1);

        tmp0 = x_ref(:,:,1:n_epo);
        tmp1 = reshape(tmp0,[size(x_ref,1),size(x_ref,2), n_epo_avg, n_obs_avg]); 
        tmp2 = permute(tmp1, [1,2,4,3]);
        tmp3 = mean(tmp2,4);
        tmp4 = permute(tmp3,[1,3,2]);
        x_ref = reshape(tmp4,size(tmp4,1),[],1);
    
        cov_ref = (x_ref * x_ref') / size(x_ref,2);
        cov_ref = cov_ref + (1e-3 * trace(cov_ref) / size(cov_ref,1) * eye(size(cov_ref,1)));    
        cov_art = (x_art * x_art') / size(x_art,2); 

        [eig_vec, eig_val] = eig(cov_art, cov_ref); % Both cov_art and cov_ref are symmetric
        eig_val = diag(eig_val);
    
        eig_val = flip(eig_val);
        eig_vec = flip(eig_vec, 2);
    
        snr_improvement = trace(eig_vec' * cov_x * eig_vec) / size(cov_x,1);
        eig_val = (eig_val - 1) / snr_improvement + 1;

        rej = eig_val > thresh;    
        if any(rej)
            pttrn = inv(eig_vec)';
            fprintf('Identified residual pulse artifacts. Rejecting %d of %d components (%.2f percent)\n', sum(rej), numel(rej), mean(rej)*100);
            R = (eye(size(cov_ref,1)) - pttrn(:,rej) * eig_vec(:, rej)');
            tmp_data = R * inEEG.data((1:inEEG.nbchan), :);
        else
            fprintf('Did not identify residual pulse artifacts\n');
            tmp_data = inEEG.data;
        end
        outEEG.data = tmp_data;
 
        errorflag = false;
    catch
        errorflag = true;
    end
end

function outcome = art_ref(inEEG,I)
    np = size(inEEG.data,1);
    I_r = size(I,1);
    I_c = size(I,2);
    tmpdata = inEEG.data;
    tmpdata = tmpdata / 1000000; % uV to V just to be consistent with MNE
    outcome = zeros(np,I_c,I_r);                
    Ivec = I(:)';
    for p = 1:np
        tmp0 = tmpdata(p,:);
        tmp1 = tmp0(Ivec);
        outcome(p,:,:) = reshape(tmp1,I_r,I_c)';
    end
end
