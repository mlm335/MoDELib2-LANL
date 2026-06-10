import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joblib, os
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

##################################################################
def get_variable_from_file(file_path, variable_name):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    for line in lines:
        stripped_line = line.strip()
        # look for lines that start with VariableName=
        if stripped_line.startswith(variable_name + "="):
            eq_index = line.find('=')
            semi_index = line.find(';', eq_index)
            if semi_index != -1:
                value_str = line[eq_index + 1:semi_index].strip()
            else:
                # fallback if no semicolon (rare but safe)
                value_str = line[eq_index + 1:].strip()
            # try to convert to number if possible
            try:
                if '.' in value_str or 'e' in value_str.lower():
                    return float(value_str)
                else:
                    return int(value_str)
            except ValueError:
                return value_str  # return as string if not numeric
    return None  # variable not found
##################################################################

eV_to_J = 1.602176634e-19
kb_eV = 8.617333262145e-5 # Boltzmann constant in eV/K


# -- RMS solute-kink interaction energy in eV -- #
def DeltaE_tilde(c):
    return 0.345*np.sqrt(c) # eV # 0.345 normal prefactor from Curtin's W-Re paper, sqrt(c) scaling from random solute distribution and central limit theorem

# -- Critical Stress Migration in Pa -- #
def tauC_tilde(c, b, a_peirls):
    return (DeltaE_tilde(c)*eV_to_J) / (b**2 * a_peirls)

# -- Migration enthalpy barrier in eV -- #
def DeltaH_km(tau, c, L, b, w, a_peirls):
    DeltaE = DeltaE_tilde(c)
    tc = tauC_tilde(c, b, a_peirls)
    term = tau / tc + 2.7 / np.sqrt(L/b)
    return DeltaE * (3.26/term + 0.035* np.sqrt(L/b) - 1.07*np.sqrt(w/b))
    # return DeltaE * (0.035* np.sqrt(L/b) - 1.07*np.sqrt(w/b))


# -- Kocks barrier -- #
def DeltaH_kocks(tau, dH0, tauP, p, q):
    x = min(tau / tauP, 1.0)          # clamp
    return dH0 * (1.0 - x**p)**q
    # return dH0 * (1 - (tau/tauP)**p)**q

def log_term(N):
    lnN = np.log(N)
    return (np.log(4*np.pi*lnN) - 4*lnN - 1.1544) / (np.sqrt(8*lnN))

def DeltaH_sol_bar(c, L, w):
    N = L/(2*w)
    return DeltaE_tilde(c) * log_term(N)

def DeltaH_dk(tau, dH0, tauP, p, q, c, L, w):
    return DeltaH_kocks(tau, dH0, tauP, p, q) + DeltaH_sol_bar(c, L, w)

# -- Activation Volume m^3 -- # (old)
def Vstar(tau_abs, c, L, b, w, a_peierls):
    tc = tauC_tilde(c, b, a_peierls)  # [Pa]
    A = 2.7 / np.sqrt(L / b)
    x = tau_abs / tc
    return 3.26 * a_peierls * b**2 * (x + A)**(-2)  # [m^3]
# =========================================================================== #

def deltaE_migration(L, b, w, c): # Zero-stress migration enthalpy ΔE in eV.
    dEtilde_p_eV = DeltaE_tilde(c)
    return dEtilde_p_eV * (1.242 * np.sqrt(L / b) - 1.07 * np.sqrt(w / b))

def Vstar_migration(L, b, a_peierls): # Zero-stress activation volume V* in m^3.
    return 0.447 * a_peierls * b * L


class MobilitySolver():
        def velocityPy(self,stress,burgers,tangent,normal,temp,dL,dt):

            # ========== Some Constants for W-Re Mobility Model ========== #
            matFilePath = "../../Library/Materials/WRe.txt" # material file path
            c_sol = get_variable_from_file(matFilePath, 'c_sol')        # [-] Re concentration
            
            # dH0_eV = get_variable_from_file(matFilePath, 'dH0_eV')      # [eV] enthalpy barrier for kink nucleation [eV]
            # p = get_variable_from_file(matFilePath, 'p')                # [-] mobility exponent
            # q = get_variable_from_file(matFilePath, 'q')                # [-] mobility exponent
            # tauP_SI = get_variable_from_file(matFilePath, 'tauC_SI')    # [Pa] Peierls stress    

            dH0_eV=1.5; 		# [eV]		enthalpy barrier for kink nucleation
            p=0.49;				# [-]		mobility exponent
            q=1.69;				# [-]		mobility exponent
            tauP_SI=2.1e9;		# [Pa]		Peierls stress (2.03e9)

            a0 = 3.1652e-10             # [m] lattice parameter in meters
            b = a0/2 * np.sqrt(3)       # [m] magnitude of Burgers vector for
            L_n = 1e-6                  # [m] length scale for nucleation
            L_m = 3e-8                  # [m] length scale for migration
            w = 12*b                    # [m] kink width in meters, estimated from literature
            X_n = (L_n - w)/b           # [m] length available for kink nucleation in meters
            X_m = L_m - w               # [m] length available for kink migration in meters
            a_peierls = 0.943*b
            a_hop = a_peierls
            h = 2.0*np.sqrt(2.0)/3.0*b   # [m] kink height, estimated from literature
            nu_D_n = 1e10               # [1/s] attempt frequency for kink nucleation in Hz, estimated from literature
            nu_D_k = 1e9              # [1/s] attempt frequency for kink migration in Hz, estimated from literature

            s_hat = burgers / np.linalg.norm(burgers)
            n_hat = normal  / np.linalg.norm(normal)
            tau = float(s_hat @ (stress @ n_hat)) * 1e9  # [Pa]
            tau_abs = abs(tau)

            # ---------------------------
            # Curtin-style barriers 
            # ---------------------------
            Hkm_eV = DeltaH_km(tau_abs, c_sol, L_m, b, w, a_peierls)
            Hdk_eV = DeltaH_dk(tau_abs, dH0_eV, tauP_SI, p, q, c_sol, L_n, w)
            Hkm_eV = max(Hkm_eV, 0.0)
            Hdk_eV = max(Hdk_eV, 0.0)
            kBT_eV = kb_eV * temp

            # ---------------------------
            # Pick Contolling Mechanism only
            # ---------------------------



            # ---------------------------
            # Dorn–Rajnak nucleation rate per length [1/(m*s)]
            # ---------------------------
            J_kp = nu_D_n * X_n * np.exp(-Hdk_eV / kBT_eV)
            t_n = 1.0 / (J_kp + 1e-300)   # [s]


            if temp == 300:
                kocks = DeltaH_kocks(tau_abs, dH0_eV, tauP_SI, p, q)
                sol_bar =  DeltaH_sol_bar(c_sol, L_n, w)
                print("T", temp, "c_sol", c_sol, "tau(GPa)", tau_abs/1e9, "Hdk_eV", Hdk_eV, "t_n", t_n)

            # ---------------------------
            # Kink drift velocity v_k (Po footnote-5 structure, no entropy) [m/s]
            # Use Joules in the sinh argument
            # ---------------------------
            kBT_J = kBT_eV * eV_to_J


            # == Original == 
            # chi = (tau * b * h * a_peierls) / (2.0 * kBT_J)  # dimensionless
            # sinh_chi = np.sinh(chi)
            # if abs(chi) < 1e-6:
            #     sinh_chi = chi # if chi is tiny, sinh(chi)~chi avoids underflow/precision loss


            # == Second Try == 
            # v_k = a_hop * nu_D_k * np.exp(-Hkm_eV / kBT_eV) # * sinh_chi  # [m/s]
            # v_k_abs = abs(v_k)
            # t_k = X_m / (2.0 * v_k_abs + 1e-300)   # [s]
            # if temp == 300:
            #     print("T", temp,  "c_sol", c_sol, "tau(GPa)", tau_abs/1e9, "Hkm_eV", Hkm_eV, "t_k", t_k)

      
            # # # == Linearization == 
            Em_eV = deltaE_migration(L_m, b, w, c_sol)
            Vstar_m3 = Vstar_migration(L_m, b, a_peierls)
            chi = (tau_abs * Vstar_m3) / (2.0 * kBT_J)        
            if chi > 50:
                sinh_chi = 0.5 * np.exp(chi)   # asymptotic
            else:
                sinh_chi = np.sinh(chi)
            v_k = 2.0 * a_hop * nu_D_k * np.exp(-Em_eV / kBT_eV) * sinh_chi  # [m/s]
            v_k_abs = abs(v_k)
            t_k = X_m / (2.0 * v_k_abs + 1e-300)   # [s]

            if temp < 200:
                print("T", temp,  "c_sol", c_sol, "tau(GPa)", tau_abs/1e9, "chi", chi, "Em(eV)", Em_eV, "t_k", t_k)


            # == Positive and Negative Drift ==
            # H_plus = DeltaH_km(tau_abs, c_sol, L_m, b, w, a_peierls)
            # H_minus = DeltaH_km(-tau_abs, c_sol, L_m, b, w, a_peierls)
            # f_plus  = nu_D_k * np.exp(-H_plus / kBT_eV)
            # f_minus = nu_D_k * np.exp(-H_minus / kBT_eV)
            # f_minus = 0
            # v_k = a_hop * (f_plus - f_minus)  # signed drift
            # v_k_abs = abs(v_k)
            # t_k = X_m / (2.0 * v_k_abs + 1e-300)
            

            # ---------------------------
            # Dislocation velocity [m/s]
            # ---------------------------
            # v = h / (t_n + t_k)

            def is_valid(x):
                return x is not None and np.isfinite(x) and x > 0

            if is_valid(t_n) and is_valid(t_k):
                # Both finite , smooth Curtin behavior
                controlling_time = t_n + t_k
            else:
                # Fallback , dominant mechanism
                valid_times = []
                if is_valid(t_n):
                    valid_times.append(t_n)
                if is_valid(t_k):
                    valid_times.append(t_k)
                if len(valid_times) > 0:
                    controlling_time = max(valid_times)
                else:
                    controlling_time = np.inf  # both invalid , frozen

            v = h / controlling_time if np.isfinite(controlling_time) else 0.0

            # print("temperature [K] =", temp)
            # print("tau [MPa] =", tau_abs/1e6)
            # print("c_sol =", c_sol)
            # print("t_n/t_k =", t_n/(t_k+1e-300))


            # print("temperature [K] =", temp)
            # print("tau [MPa] =", tau_abs/1e6)
            # print("Hdk [eV] =", Hdk_eV, "Hsol [eV] =", DeltaH_sol_bar(c_sol,L,w), "H0 [eV] =", DeltaH_kocks(tau_abs,dH0_eV,tauP_SI,p,q))
            # print("Hkm [eV] =", Hkm_eV)
            # print("kBT [eV] =", kb_eV*temp)
            # print("exp(-Hdk/kBT) =", np.exp(-Hdk_eV/(kb_eV*temp)))
            # print("exp(-Hkm/kBT) =", np.exp(-Hkm_eV/(kb_eV*temp)))
            # print("J_kp*X =", (nu_D/w)*np.exp(-Hdk_eV/(kb_eV*temp))*X_n)
            # print("t_n =", t_n, "t_k =", t_k)
            # print("chi =", chi, "sinh(chi) =", np.sinh(chi))
            # print("v_k =", v_k)

            return v

