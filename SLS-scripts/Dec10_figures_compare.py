from Nov21_Cart_synthesis_Qparam import (
    plantMdelta, Qsyn, Q_synthesis, Ts)
import SLSsyn as Ss
import control as co
import matplotlib.pyplot as plt
import numpy as np
import pickle

# 
COLOR_idl = 'C0'
COLOR_ad = 'C1'
COLOR_Q = 'C2'
COLOR_bnd = 'gray'


def load_Qparam_controller(name):
    with open(name, 'rb') as f:
        K = pickle.load(f)
    return K


def sim_step(sys, T_sim):
    _, y = co.step_response(sys, T_sim)
    return y[0, :]


def fig_admittance_shaping():
    Pz_nom = plantMdelta(
        E_gain=50, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)
    Pz_contracted = plantMdelta(
        E_gain=100, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)
    Pz_relaxed = plantMdelta(
        E_gain=20, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=2, N_out=2)

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
    fig = plt.figure(figsize=[7, 3.5])
    ax0 = fig.add_subplot(2, 3, 1)
    ax1 = fig.add_subplot(2, 3, 2, sharex=ax0)
    ax2 = fig.add_subplot(2, 3, 3, sharex=ax0)

    # time response
    y0 = sim_step(H_relaxed_ideal, T_sim)
    y1 = sim_step(H_nom_ideal, T_sim)
    y2 = sim_step(H_contracted_ideal, T_sim)
    ax0.plot(T_sim, y0, label='idl.', linestyle='-', color=COLOR_idl)
    ax1.plot(T_sim, y1, label='idl.', linestyle='-', color=COLOR_idl)
    ax2.plot(T_sim, y2, label='idl.', linestyle='-', color=COLOR_idl)

    y0 = sim_step(H_relaxed_ad[0, 2], T_sim)
    y1 = sim_step(H_nom_ad[0, 2], T_sim)
    y2 = sim_step(H_contracted_ad[0, 2], T_sim)
    ax0.plot(T_sim, y0, label='ad', linestyle='--', color=COLOR_ad)
    ax1.plot(T_sim, y1, label='ad', linestyle='--', color=COLOR_ad)
    ax2.plot(T_sim, y2, label='ad', linestyle='--', color=COLOR_ad)

    y0 = sim_step(H_relaxed_Q[0, 2], T_sim)
    y1 = sim_step(H_nom_Q[0, 2], T_sim)
    y2 = sim_step(H_contracted_Q[0, 2], T_sim)
    ax0.plot(T_sim, y0, label='Q', linestyle='-.', color=COLOR_Q)
    ax1.plot(T_sim, y1, label='Q', linestyle='-.', color=COLOR_Q)
    ax2.plot(T_sim, y2, label='Q', linestyle='-.', color=COLOR_Q)

    # # ax0.legend(ncol=3, loc='best', framealpha=1)
    # for ax in [ax0, ax1, ax2]:
    #     ax.set_xlabel('Time (sec)')
    # ax0.set_ylabel('Velocity (m/s)')

    # freq-response
    freqs = np.logspace(-2, np.log10(np.pi / Ts), 100)

    ax3 = fig.add_subplot(2, 3, 4)
    ax4 = fig.add_subplot(2, 3, 5, sharex=ax3, sharey=ax3)
    ax5 = fig.add_subplot(2, 3, 6, sharex=ax3, sharey=ax3)
    for ax in [ax3, ax4, ax5]:
        ax.set_xscale('log')
        ax.set_xlim([3e-1, 60])

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
 
    # hide x labels
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)
   
    plt.tight_layout()
    plt.savefig("Dec10_admittance_shaping.pdf")
    plt.show()


def fig_noise_attenuation():
    Pz_nom = plantMdelta(
        E_gain=50, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)
    Pz_contracted = plantMdelta(
        E_gain=100, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)
    Pz_relaxed = plantMdelta(
        E_gain=20, wI=1.0, sys_type='33_mult_unt', m_int=0.1, N_in=1, N_out=1)

    noise_atten_func = Qsyn.lambda_log_interpolate(
        [[0.1, 0.1], [25, 0.1], [25, 0.008], [100, 0.008]])

    K_ad = co.c2d(co.tf([1], [2.5, 12, 0]), Ts)
    K_Q = load_Qparam_controller('Nov21_ctrl.pkl')

    # pre-processing
    H_nom_ad = Ss.lft(Pz_nom, K_ad)
    H_relaxed_ad = Ss.lft(Pz_relaxed, K_ad)
    H_contracted_ad = Ss.lft(Pz_contracted, K_ad)

    H_nom_Q = Ss.lft(Pz_nom, K_Q)
    H_relaxed_Q = Ss.lft(Pz_relaxed, K_Q)
    H_contracted_Q = Ss.lft(Pz_contracted, K_Q)

    # fig init
    fig = plt.figure(figsize=[7, 2.5])
    ax0 = fig.add_subplot(2, 3, 1)
    ax1 = fig.add_subplot(2, 3, 2, sharex=ax0, sharey=ax0)
    ax2 = fig.add_subplot(2, 3, 3, sharex=ax0, sharey=ax0)
    ax3 = fig.add_subplot(2, 3, 4, sharex=ax0)
    ax4 = fig.add_subplot(2, 3, 5, sharex=ax0, sharey=ax3)
    ax5 = fig.add_subplot(2, 3, 6, sharex=ax0, sharey=ax3)

    # hide x labels
    plt.setp(ax0.get_xticklabels(), visible=False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax1.get_yticklabels(), visible=False)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    plt.setp(ax5.get_yticklabels(), visible=False)

    ax3.set_yticks(np.pi * np.r_[-0.5, 0, 0.5])
    ax3.set_yticklabels(['$-\pi/2$', '$0$', '$\pi/2$'])

    # freq-response
    freqs = np.linspace(1e-3, np.pi / Ts, 1000)
    for ax in [ax0, ax1, ax2]:
        # ax.set_xscale('log')
        ax.set_yscale('log')

    for ax in [ax3, ax4, ax5]:
        # ax.plot([freqs[0], freqs[-1]], [np.pi / 2, np.pi / 2], '-', c='gray')
        # ax.plot([freqs[0], freqs[-1]], [- np.pi / 2, - np.pi / 2], '-', c='gray')
        ax.grid()
        ax.set_xticks([0, 10, 20, 30, 40])

    mag0, phase0, _ = H_relaxed_ad[0, 0].freqresp(freqs)
    mag1, phase1, _ = H_nom_ad[0, 0].freqresp(freqs)
    mag2, phase2, _ = H_contracted_ad[0, 0].freqresp(freqs)
    ax0.plot(freqs, mag0[0, 0], label='ad.', linestyle='--', c=COLOR_ad)
    ax1.plot(freqs, mag1[0, 0], label='ad.', linestyle='--', c=COLOR_ad)
    ax2.plot(freqs, mag2[0, 0], label='ad.', linestyle='--', c=COLOR_ad)
    ax3.plot(freqs, phase0[0, 0], label='ad.', linestyle='--', c=COLOR_ad)
    ax4.plot(freqs, phase1[0, 0], label='ad.', linestyle='--', c=COLOR_ad)
    ax5.plot(freqs, phase2[0, 0], label='ad.', linestyle='--', c=COLOR_ad)

    mag0, phase0, _ = H_relaxed_Q[0, 0].freqresp(freqs)
    mag1, phase1, _ = H_nom_Q[0, 0].freqresp(freqs)
    mag2, phase2, _ = H_contracted_Q[0, 0].freqresp(freqs)
    ax0.plot(freqs, mag0[0, 0], label='Q', linestyle='-.', c=COLOR_Q)
    ax1.plot(freqs, mag1[0, 0], label='Q', linestyle='-.', c=COLOR_Q)
    ax2.plot(freqs, mag2[0, 0], label='Q', linestyle='-.', c=COLOR_Q)
    ax3.plot(freqs, phase0[0, 0], label='Q', linestyle='-.', c=COLOR_Q)
    ax4.plot(freqs, phase1[0, 0], label='Q', linestyle='-.', c=COLOR_Q)
    ax5.plot(freqs, phase2[0, 0], label='Q', linestyle='-.', c=COLOR_Q)

    mag3 = noise_atten_func(freqs)
    ax0.plot(freqs, mag3, label='bnd.', linestyle='-', c=COLOR_bnd)
    ax1.plot(freqs, mag3, label='bnd.', linestyle='-', c=COLOR_bnd)
    ax2.plot(freqs, mag3, label='bnd.', linestyle='-', c=COLOR_bnd)

    # post-process
    for ax in [ax0, ax1, ax2]:
        ax.grid(which='major')
        ax.set_xlim(0, 45)
        ax.set_ylim(0.005, 0.15)
    ax3.set_ylim(-2.2, 2.2)
    plt.tight_layout()
    plt.savefig("Dec10_noise.pdf")
    plt.show()


if __name__ == '__main__':
    fig_admittance_shaping()
    # fig_noise_attenuation()
