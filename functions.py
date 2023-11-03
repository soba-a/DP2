# ALL FUNCTIONS IN METRIC UNITS!!!!!!!!!!!!!!
import numpy as np

def calc_Se(Se_prime, ka, kb, kc, kd, ke, kf): # kf = 1 for this class!
    return Se_prime * ka * kb * kc * kd * ke * kf

def calc_Kb(d):
    if d < 7.62: return 1
    elif 7.62 <= d <= 51: return 1.24 * d**(-0.107)
    elif 51 < d < 254: return 1.51 * d**(-0.157)
    else: return -1

def calc_K_f(K_t, d, S_ut):
    if 340 <= S_ut <= 1700: #MPa
        a = 1.24 - 2.25*10**(-3)*S_ut + 1.6*10**(-6)*S_ut**2 - 4.1*10**(-10)*S_ut**3
        return 1 + (K_t - 1)/(1 + (a/(d/2))**(0.5))
    else:
        return -1

def calc_K_fs(K_ts, d, S_ut):
    if 340 <= S_ut <= 1700: # MPa
        a = 0.958 - 1.83*10**(-3)*S_ut + 1.43*10**(-6)*S_ut**2 - 4.11*10**(-10)*S_ut**3
        return 1 + (K_ts - 1)/(1 + (a/(d/2))**(0.5))
    else:
        return -1
    
def sig_a_prime(K_f, M, d):
    return K_f * 32 * M / (np.pi * d**3)

def sig_m_prime(K_fs, T, d):
    return 3**(-0.5) * K_fs * 16 * T / (np.pi * d**3)

def goodman_safety(sig_a, sig_m, S_e, S_ut):
    return (sig_a / S_e + sig_m / S_ut)**(-1)

def asme_safety(sig_a, sig_m, S_e, S_y):
    return (sig_a / S_e + sig_m / S_y)**(-0.5)

def yield_check(sig_a, sig_m, S_y):
    return S_y / (sig_a + sig_m)