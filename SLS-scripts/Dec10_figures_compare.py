from Nov21_Cart_synthesis_Qparam import (
    plantMdelta, Qsyn, Q_synthesis, Ts)
import SLSsyn as Ss
import control as co
import matplotlib.pyplot as plt
import numpy as np
import pickle


def load_Qparam_controller(name):
    with open(name, 'rb') as f:
        K = pickle.load(f)
    return K


def sim_step(sys, T_sim):
    _, y = co.step_response(sys, T_sim)
    return y[0, :]


def fig_admittance_shaping():
    Pz_nom = plantMdelta(
        E_gain=50, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)
    Pz_contracted = plantMdelta(
        E_gain=100, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)
    Pz_relaxed = plantMdelta(
        E_gain=20, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)

    m_des = 2.5
    b_des = 12
    K_ad = co.c2d(co.tf([1], [2.5, 12, 0]), Ts)
    K_Q = load_Qparam_controller('Nov21_ctrl.pkl')

    T_sim = np.arange(600) * Ts

    # pre-processing
    H_nom_ad = Ss.lft(Pz_nom, K_ad)
    H_relaxed_ad = Ss.lft(Pz_relaxed, K_ad)
    H_contracted_ad = Ss.lft(Pz_contracted, K_ad)

    H_nom_Q = Ss.lft(Pz_nom, K_Q)
    H_relaxed_Q = Ss.lft(Pz_relaxed, K_Q)
    H_contracted_Q = Ss.lft(Pz_contracted, K_Q)

    H_nom_ideal = co.c2d(co.tf([50, 0], [m_des, b_des, 0 + 50]), Ts)
    H_relaxed_ideal = co.c2d(co.tf([20, 0], [m_des, b_des, 0 + 20]), Ts)
    H_contracted_ideal = co.c2d(co.tf([100, 0], [m_des, b_des, 0 + 100]), Ts)

    # fig init
    fig = plt.figure(figsize=[8, 4])
    ax0 = fig.add_subplot(2, 3, 1)
    ax1 = fig.add_subplot(2, 3, 2, sharex=ax0)
    ax2 = fig.add_subplot(2, 3, 3, sharex=ax0)

    # time response
    y0 = sim_step(H_relaxed_ideal, T_sim)
    y1 = sim_step(H_nom_ideal, T_sim)
    y2 = sim_step(H_contracted_ideal, T_sim)
    ax0.plot(T_sim, y0, label='idl.', linestyle='-')
    ax1.plot(T_sim, y1, label='idl.', linestyle='-')
    ax2.plot(T_sim, y2, label='idl.', linestyle='-')

    y0 = sim_step(H_relaxed_ad[0, 2], T_sim)
    y1 = sim_step(H_nom_ad[0, 2], T_sim)
    y2 = sim_step(H_contracted_ad[0, 2], T_sim)
    ax0.plot(T_sim, y0, label='ad', linestyle='--')
    ax1.plot(T_sim, y1, label='ad', linestyle='--')
    ax2.plot(T_sim, y2, label='ad', linestyle='--')

    y0 = sim_step(H_relaxed_Q[0, 2], T_sim)
    y1 = sim_step(H_nom_Q[0, 2], T_sim)
    y2 = sim_step(H_contracted_Q[0, 2], T_sim)
    ax0.plot(T_sim, y0, label='Q', linestyle='-.')
    ax1.plot(T_sim, y1, label='Q', linestyle='-.')
    ax2.plot(T_sim, y2, label='Q', linestyle='-.')

    # ax0.legend(ncol=3, loc='best', framealpha=1)
    for ax in [ax0, ax1, ax2]:
        ax.set_xlabel('Time (sec)')
    ax0.set_ylabel('Velocity (m/s)')

    # freq-response
    freqs = np.logspace(-2, np.log10(np.pi / Ts), 100)

    ax3 = fig.add_subplot(2, 3, 4)
    ax4 = fig.add_subplot(2, 3, 5, sharex=ax3, sharey=ax3)
    ax5 = fig.add_subplot(2, 3, 6, sharex=ax3, sharey=ax3)
    for ax in [ax3, ax4, ax5]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim([1e-3, 30])
        ax.set_xlabel('Frequency (rad/s)')
    ax3.set_ylabel('Magnitude')

    mag0, _, _ = H_relaxed_ideal.freqresp(freqs)
    mag1, _, _ = H_nom_ideal.freqresp(freqs)
    mag2, _, _ = H_contracted_ideal.freqresp(freqs)
    ax3.plot(freqs, mag0[0, 0], label='idl.', linestyle='-')
    ax4.plot(freqs, mag1[0, 0], label='idl.', linestyle='-')
    ax5.plot(freqs, mag2[0, 0], label='idl.', linestyle='-')

    mag0, _, _ = H_relaxed_ad[0, 2].freqresp(freqs)
    mag1, _, _ = H_nom_ad[0, 2].freqresp(freqs)
    mag2, _, _ = H_contracted_ad[0, 2].freqresp(freqs)
    ax3.plot(freqs, mag0[0, 0], label='ad.', linestyle='--')
    ax4.plot(freqs, mag1[0, 0], label='ad.', linestyle='--')
    ax5.plot(freqs, mag2[0, 0], label='ad.', linestyle='--')

    mag0, _, _ = H_relaxed_Q[0, 2].freqresp(freqs)
    mag1, _, _ = H_nom_Q[0, 2].freqresp(freqs)
    mag2, _, _ = H_contracted_Q[0, 2].freqresp(freqs)
    ax3.plot(freqs, mag0[0, 0], label='Q', linestyle='-.')
    ax4.plot(freqs, mag1[0, 0], label='Q', linestyle='-.')
    ax5.plot(freqs, mag2[0, 0], label='Q', linestyle='-.')

    # post-process
    for ax in [ax0, ax1, ax2, ax3, ax4, ax5]:
        ax.grid(which='major')
    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    fig_admittance_shaping()
