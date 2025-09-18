import numpy as np
import pandas as pd
import joblib

## Load model
loaded_model = joblib.load('/Users/matthewmaron/Documents/MoDELib2/python/MobilityLaw/GBR_model.pkl')

def sort_eigenvalues_and_vectors(eigenvalues, eigenvectors):
    idx = eigenvalues.argsort()[::-1]  # Sort in descending order
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    return eigenvalues, eigenvectors


class MobilitySolver():
        def __init__(self, log_path="MobilityLog.xlsx", log_sheet="Sheet1"):
            self.log_path = log_path
            self.log_sheet = log_sheet

        def velocityPy(self,stress,burgers,tangent,normal,temp,dL,dt):

            # --- Missorientation angle in degrees ---
            b_hat = burgers / np.linalg.norm(burgers)
            t_hat = tangent / np.linalg.norm(tangent)
            cos_theta = np.clip(np.dot(b_hat, t_hat), -1.0, 1.0)
            miss_angle = float(np.degrees(np.arccos(cos_theta)))

            # --- Eigenvalues and eigenvectors, Project eigenvectors along tangent, Compute absolute values ---
            eigValue, eigVector = np.linalg.eig(stress)
            eigValue, eigVector = sort_eigenvalues_and_vectors(eigValue, eigVector)
            projection = eigVector.T @ tangent
            avg_up_for_sim = pd.DataFrame(index=range(1), columns=['nev1', 'nev2', 'nev3'])
            avg_up_for_sim.loc[:, 'nev1'] = projection[0]
            avg_up_for_sim.loc[:, 'nev2'] = projection[1]
            avg_up_for_sim.loc[:, 'nev3'] = projection[2]
            avg_up_for_sim['anev1'] = avg_up_for_sim['nev1'].abs()
            avg_up_for_sim['anev2'] = avg_up_for_sim['nev2'].abs()
            avg_up_for_sim['anev3'] = avg_up_for_sim['nev3'].abs()

            # --- Feature Set ---
            features = pd.DataFrame({
                'e1': eigValue[0],
                'e2': eigValue[1],
                'e3': eigValue[2],
                'anev1': avg_up_for_sim['anev1'],
                'anev2': avg_up_for_sim['anev2'],
                'anev3': avg_up_for_sim['anev3'],
                'temp': temp,
                'miss': miss_angle
            })

            # --- predict velocity (Å/ps -> m/s) ---
            velocity_pred = loaded_model.predict(features)
            velocity_mps = float(velocity_pred[0]) * 1000.0

            return velocity_mps
