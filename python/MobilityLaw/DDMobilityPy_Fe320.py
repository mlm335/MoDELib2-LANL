import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import joblib, os
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

## Load model
loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2/python/MobilityLaw/gpr_8_feature_sqrt_RBF_1001.joblib')

def sort_eigenvalues_and_vectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]  # Sort in descending order
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors

def angle_atan2(b_unit, n_unit, t_unit, f_unit):
    # Compute the angle using atan2 for better numerical stability
    cos_theta = np.sum(b_unit*t_unit, axis=-1)
    sin_theta = np.sum(np.cross(b_unit, t_unit)*n_unit, axis=-1)
    th_deg = np.degrees(np.arctan2(sin_theta, cos_theta))
    return (th_deg % 180 + 180) % 180  # Ensure angle is in [0, 180)

class MobilitySolver():
        def velocityPy(self,stress,burgers,tangent,normal,temp,dL,dt):
                #Missorientation angle in degrees
                # miss_product = np.dot(burgers, tangent)
                # magnitude_a = np.linalg.norm(burgers)
                # magnitude_b = np.linalg.norm(tangent)
                # cos_theta = np.clip(miss_product / (magnitude_a * magnitude_b), -1, 1)
                # angle_radians = np.arccos(cos_theta)
                # miss_angle = np.degrees(angle_radians)


                f_pk = np.cross(stress@burgers, tangent)
                f_g = f_pk - np.dot(f_pk,normal)*normal
                f_g_norm = f_g / np.linalg.norm(f_g)
                miss_angle = angle_atan2(burgers, normal, tangent, f_g_norm)
                
                #Eigenvalues and eigenvectors, Project eigenvectors along tangent, Compute absolute values
                eigValue, eigVector = np.linalg.eig(stress)
                eigValue, eigVector = sort_eigenvalues_and_vectors(eigValue, eigVector)
                projection = eigVector.T @ tangent
                anev = np.abs(projection)

                # avg_up_for_sim = pd.DataFrame(index=range(1), columns=['nev1', 'nev2', 'nev3'])
                # avg_up_for_sim.loc[:, 'nev1'] = projection[0]
                # avg_up_for_sim.loc[:, 'nev2'] = projection[1]
                # avg_up_for_sim.loc[:, 'nev3'] = projection[2]
                # avg_up_for_sim['anev1'] = avg_up_for_sim['nev1'].abs()
                # avg_up_for_sim['anev2'] = avg_up_for_sim['nev2'].abs()
                # avg_up_for_sim['anev3'] = avg_up_for_sim['nev3'].abs()

                # Feature Set
                features = pd.DataFrame({
                    'e1': eigValue[0],
                    'e2': eigValue[1],
                    'e3': eigValue[2],
                    'anev1': anev[0],
                    'anev2': anev[1],
                    'anev3': anev[2],
                    'temp': temp,
                    'miss': miss_angle
                }, index=[0])

                # features = pd.DataFrame([{
                #         'e1': eigValue[0],
                #         'e2': eigValue[1],
                #         'e3': eigValue[2],
                #         'temp': temp,
                #         'miss': miss_angle
                #         }])
                
                velocity_pred = loaded_model.predict(features)
                velocity_mps = float(velocity_pred) 
                return velocity_mps**2 
