from typing import Tuple

import numpy as np

from scipy.linalg import eigh
import matplotlib.pyplot as plt
from pyriemann.utils.mean import mean_riemann
import mne

from mne.decoding.base import BaseEstimator
from mne.decoding.mixin import TransformerMixin
import warnings

class EventBasedArtifactProjector(TransformerMixin, BaseEstimator):

    def __init__(
        self,
        event_id : str,
        threshold : float = 1.15,
        events_per_erp : int = 25,
        n_repetitions : int = 1,
        offset : float = -0.2,
        duration : float = 0.7,
        picks : str = 'eeg',
    ) -> None:
        super().__init__() 
        self._event_id = event_id  
        self._threshold = threshold
        self._events_per_erp = events_per_erp
        self._n_repetitions = n_repetitions
        self._offset = offset
        self._duration = duration
        self._info = None
        self._picks = picks

    def fit(self, raw : mne.io.Raw, y = None):

        ecg_evt_ixs = np.flatnonzero(raw.annotations.description == self._event_id)
        if ecg_evt_ixs.shape[0] == 0:
            warnings.warn(f'Did not find event "{self._event_id}" in raw object.')
            return self

        ecg_trig = np.round((raw.annotations[ecg_evt_ixs].onset -raw.first_time) * raw.info['sfreq']).astype(int)    

        offset = int(self._offset*raw.info['sfreq'])
        duration = int(self._duration*raw.info['sfreq'])

        self._info = raw.copy().pick(self._picks).info

        # create slicing mask
        I = ecg_trig[:,None] + np.arange(duration)[None,:] + offset
        # remove not feasible rows
        I = I[((I >= 0) & (I < raw.times.shape[-1])).all(axis=1),:]

        # get artifact data
        x_art = raw.get_data(picks=self._picks)[:,I].transpose((0,2,1))

        n_chan = x_art.shape[0]
        n_events = x_art.shape[-1]

        if n_events < self._events_per_erp:
            warnings.warn(f'Too few events (n={n_events}) in the data. Consider fewer events to compute ERPs ("adjust events_per_erp parameter.").')
            return self

        # compute beat ERPs and the spatial covariance matrices per ERP
        n_erp = x_art.shape[2] // self._events_per_erp
        n_events = n_erp * self._events_per_erp
        x_art = x_art[:,:,:n_events].reshape((n_chan, duration, n_erp, self._events_per_erp)).mean(axis=3).transpose((2, 0, 1))

        cov_art = (x_art @ x_art.swapaxes(-2,-1)) / x_art.shape[-1]
        cov_shrinkage = 1e-3 * np.median(np.trace(cov_art, axis1=-2,axis2=-1)) / cov_art.shape[-1] * np.eye(cov_art.shape[-1])
        cov_art += cov_shrinkage
        # average all available ERPs
        cov_art = mean_riemann(cov_art)

        covs_ref = []
        covs_x = []

        for _ in range(self._n_repetitions):

            rnd_bins = np.random.randint(0, raw.times.shape[-1] - duration, ecg_trig.shape[0])
            IR = rnd_bins[:,None] + np.arange(duration)[None,:]
            # remove not feasible rows
            IR = IR[((IR >= 0) & (IR < raw.times.shape[-1])).all(axis=1),:]
            IR = IR[:I.shape[0]]


            x_ref = raw.get_data(picks=self._picks)[:,IR].transpose((0,2,1))
            x_ref = x_ref[:,:,:n_events].reshape((n_chan, duration, n_erp, self._events_per_erp)).transpose((2, 3, 0, 1))
            cov_x = (x_ref @ x_ref.swapaxes(-2,-1)) / x_ref.shape[-1]
            cov_x += cov_shrinkage
            cov_x = mean_riemann(cov_x.mean(axis=1))
            x_ref_avg = x_ref.mean(axis=1)

            cov_ref = (x_ref_avg @ x_ref_avg.swapaxes(-2,-1)) / x_ref_avg.shape[-1]
            cov_ref += cov_shrinkage 
            cov_ref = mean_riemann(cov_ref)

            covs_ref.append(cov_ref)
            covs_x.append(cov_x)

        covs_ref = np.stack(covs_ref)
        covs_x = np.stack(covs_x)

        cov_ref = mean_riemann(covs_ref)
        cov_x = mean_riemann(covs_x)

        eig_val, eig_vec = eigh(cov_art, cov_ref)
        # sort in descending order
        eig_val = eig_val[::-1]
        eig_vec = eig_vec[:,::-1]
        self._snr_improvement = np.trace(eig_vec.T @ cov_x @ eig_vec)/cov_x.shape[0]
        self._eig_vec = eig_vec
        self._eig_val = (eig_val - 1) / self._snr_improvement + 1
        self._rej = self._eig_val > self._threshold
        self._pttrn  = np.linalg.inv(eig_vec).T
        print(f"Identified artifacts. Rejecting {self._rej.sum()} of {len(self._rej)} components ({self._rej.mean()*100:.0f}%).")

        return self

    def transform(self, raw : mne.io.Raw, y = None):

        chan_idxs = mne.channel_indices_by_type(raw.info)[self._picks]
        raw_ch_names = [raw.ch_names[ix] for ix in chan_idxs]
        assert len(raw_ch_names) == len(self._info['ch_names'])
        assert all([trans_name == fit_name for trans_name, fit_name in zip(raw_ch_names, self._info['ch_names'])])

        n_chan = len(chan_idxs)
        R = (np.eye(n_chan) - self._pttrn[:,self._rej] @ self._eig_vec[:, self._rej].T)

        raw._data[chan_idxs,:] = R @ raw._data[chan_idxs, :]
        return raw

    def plot(self, n_components = None):

        if n_components is None:
            n_components = int(self._rej.sum())
        if n_components > 0:
            fig, axes = plt.subplots(1, n_components, figsize=(n_components*2,2))
            axes = axes if n_components > 1 else [axes]
            for ix, ax in enumerate(axes):
                _ = mne.viz.plot_topomap(self._pttrn[:,ix], self._info, axes=ax, show=False)
                ax.set_title(f'eigval={self._eig_val[ix]:.2f}\nsnr_impr.={self._snr_improvement:.2f}\nrej={self._rej[ix]}')
            return fig
        else:
            return None
