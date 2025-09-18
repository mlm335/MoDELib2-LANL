import numpy as np
from scipy.interpolate import griddata

def miller_to_latex(vec):
    def fmt(i):
        return f"\\bar{{{abs(i)}}}" if i < 0 else str(i)
    h, k, l = map(int, vec)
    return f"$[{fmt(h)}{fmt(k)}{fmt(l)}]$"


def angleAxis(theta,a):
    a=a/np.linalg.norm(a) # normalize axis a
    return np.identity(3)*np.cos(theta)+np.array([[0,-a[2],a[1]],[a[2],0,-a[0]],[-a[1],a[0],0]])*np.sin(theta)+np.outer(a,a)*(1-np.cos(theta))


def normalize(v):
    v = np.array(v, dtype=float)
    norm = np.linalg.norm(v)
    if norm == 0:
        raise ValueError("Zero vector cannot be normalized.")
    return v / norm


def standard_directions():
    directions = [
        [1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],
        [0,1,1],[0,-1,1],[-1,0,1],[1,0,1],[-1,-1,0],[1,-1,0],[-1,1,0],[1,1,0],
        [1,1,1],[-1,1,1],[-1,-1,1],[1,-1,1]]
    return [(normalize(np.array(direction)),np.array(direction)) for direction in directions]

def midpoint_directions():
    directions = [
        [-2, -2, 1], [-2, -1, 0], [-2, -1, 1], [-2, -1, 2], [-2, 0, 1],
        [-2, 1, 0], [-2, 1, 1], [-2, 1, 2], [-2, 2, 1], [-1, -2, 0],
        [-1, -2, 1], [-1, -2, 2], [-1, -1, 2], [-1, 0, 2], [-1, 1, 2],
        [-1, 2, 0], [-1, 2, 1], [-1, 2, 2], [0, -2, 1], [0, -1, 2],
        [0, 1, 2], [0, 2, 1], [1, -2, 0], [1, -2, 1], [1, -2, 2],
        [1, -1, 2], [1, 0, 2], [1, 1, 2], [1, 2, 0], [1, 2, 1],
        [1, 2, 2], [2, -2, 1], [2, -1, 0], [2, -1, 1], [2, -1, 2],
        [2, 0, 1], [2, 1, 0], [2, 1, 1], [2, 1, 2], [2, 2, 1]]
    return [(normalize(np.array(direction)), np.array(direction)) for direction in directions]


def stereographic_projection(v):
    v = normalize(v)
    x, y, z = v
    return x / (1 + z), y / (1 + z)


def plot_standard_great_circles(ax, num_points=1500, color='k', linewidth=2):
    great_circle_planes = {
        "(100)": [1, 0, 0],
        "(010)": [0, 1, 0],
        "(110)": [1, 1, 0],
        "(1-10)": [1, -1, 0],
        "(101)": [1, 0, 1],
        "(011)": [0, 1, 1],
        "(01-1)": [0, -1, 1],
        "(10-1)": [-1, 0, 1]
    }
    for normal in great_circle_planes.values():
        n = normalize(np.array(normal))
        a = np.array([1, 0, 0])
        if np.allclose(np.abs(np.dot(n, a)), 1.0):
            a = np.array([0, 1, 0])
        v1 = normalize(np.cross(n, a))
        v2 = normalize(np.cross(n, v1))
        angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        points = [v1 * np.cos(t) + v2 * np.sin(t) for t in angles]
        proj_points = [stereographic_projection(p) for p in points]
        rotated = [(y, -x) for (x, y) in proj_points]
        x, y = zip(*[(x, y) for (x, y) in rotated if x**2 + y**2 <= 1])
        ax.plot(x, y, color=color, linewidth=linewidth)


def get_bcc_slip_systems():
    # Format: (slip_direction, slip_plane_normal)
    raw_slip_systems = [
        ([1, -1, 1], [1, 1, 0]),
        ([-1, 1, 1], [1, 1, 0]),
        ([1, 1, 1], [-1, 1, 0]),
        ([-1, -1, 1], [-1, 1, 0]),
        ([-1, -1, 1], [1, 0, 1]),
        ([-1, 1, 1], [1, 0, 1]),
        ([1, 1, 1], [-1, 0, 1]),
        ([1, -1, 1], [-1, 0, 1]),
        ([-1, -1, 1], [0, 1, 1]),
        ([1, -1, 1], [0, 1, 1]),
        ([1, 1, 1], [0, -1, 1]),
        ([-1, 1, 1], [0, -1, 1]) ]
    return [(normalize(np.array(s)), normalize(np.array(m))) for s, m in raw_slip_systems]


def compute_full_stereographic_velocity_map(mobility_obj, slip_systems, T=300, s0=0.01, theta=0, res=500):
    # Generate 2D stereographic grid on unit circle
    q1_vals = np.linspace(-1, 1, res)
    q2_vals = np.linspace(-1, 1, res)
    q1_grid, q2_grid = np.meshgrid(q1_vals, q2_vals)
    mask = q1_grid**2 + q2_grid**2 <= 1
    q1_flat = q1_grid[mask]
    q2_flat = q2_grid[mask]
    r2 = q1_flat**2 + q2_flat**2
    # Inverse stereographic projection to 3D
    x = 2 * q1_flat / (1 + r2)
    y = 2 * q2_flat / (1 + r2)
    z = (1 - r2) / (1 + r2)
    directions = np.stack((x, y, z), axis=1)
    # Compute max velocity per direction
    v_values = []
    for m in directions:
        m = m / np.linalg.norm(m)
        S = np.outer(m, m) * s0
        max_v = 0.0
        for b, n in slip_systems:
            b = b / np.linalg.norm(b)
            n = n / np.linalg.norm(n)
            xi=angleAxis(theta,n)@b # rotate b about n to define the line tangent
            try:
                v = mobility_obj.velocity(S, b, xi, n, T)
                if v > max_v:
                    max_v = v
            except Exception as e:
                print("Error:", e)
                continue
        v_values.append(max_v)
    v_values = np.array(v_values)
    # Interpolate to full grid
    grid_z = griddata((q1_flat, q2_flat), v_values, (q1_grid, q2_grid), method='cubic')
    grid_z[~mask] = np.nan
    return grid_z, q1_grid, q2_grid
