"""One satellite, one instant: Earth, satellite, and the satellite->Sun vector (plotly).

sat, sun : Cartesian [x, y, z] in km, Earth-centred, same frame.
           No coordinate conversion happens here; convert before calling.

Should render (via plotly) on localhost as an interactive web app
"""
#! THIS MODULE NEEDS TO BE REMOVED BEFORE MERGING
import numpy as np
import plotly.graph_objects as go

R_EARTH = 6371.0  # km


def _arrow(fig, start, vec, color, name, width=6):
    """Line from `start` along `vec` with a cone head at the tip."""
    end = start + vec
    fig.add_trace(go.Scatter3d(
        x=[start[0], end[0]], y=[start[1], end[1]], z=[start[2], end[2]],
        mode="lines", line=dict(color=color, width=width),
        name=name, legendgroup=name, hoverinfo="skip"))
    head = 0.18 * np.linalg.norm(vec)
    d = vec / np.linalg.norm(vec)
    fig.add_trace(go.Cone(
        x=[end[0]], y=[end[1]], z=[end[2]], u=[d[0]], v=[d[1]], w=[d[2]],
        anchor="tip", sizemode="absolute", sizeref=head,
        colorscale=[[0, color], [1, color]], showscale=False,
        name=name, legendgroup=name, showlegend=False, hoverinfo="skip"))


def plot_sat_sun(sat, sun, R=R_EARTH, arrow_len=None, show=True):
    sat_xyz = np.asarray(sat, dtype=float).reshape(3)
    sun_xyz = np.asarray(sun, dtype=float).reshape(3)
    sat_to_sun = sun_xyz - sat_xyz                 # the vector subtraction
    u = sat_to_sun / np.linalg.norm(sat_to_sun)    # unit direction

    r_sat = np.linalg.norm(sat_xyz)
    altitude = r_sat - R
    zenith = sat_xyz / r_sat
    angle = np.degrees(np.arccos(np.clip(u @ zenith, -1, 1)))

    fig = go.Figure()

    # Earth (opaque, lit sphere)
    a = np.linspace(0, 2 * np.pi, 120)
    b = np.linspace(0, np.pi, 60)
    fig.add_trace(go.Surface(
        x=R * np.outer(np.cos(a), np.sin(b)),
        y=R * np.outer(np.sin(a), np.sin(b)),
        z=R * np.outer(np.ones_like(a), np.cos(b)),
        colorscale=[[0, "#2b5fbf"], [1, "#2b5fbf"]], showscale=False,
        lighting=dict(ambient=0.55, diffuse=0.8, specular=0.1, roughness=0.9),
        hoverinfo="skip", name="Earth", showlegend=True))

    # Satellite, plus a stem from the sub-satellite point to show altitude
    foot = zenith * R
    fig.add_trace(go.Scatter3d(
        x=[foot[0], sat_xyz[0]], y=[foot[1], sat_xyz[1]], z=[foot[2], sat_xyz[2]],
        mode="lines", line=dict(color="white", width=4, dash="dash"),
        name=f"altitude {altitude:.1f} km", hoverinfo="skip"))
    fig.add_trace(go.Scatter3d(
        x=[sat_xyz[0]], y=[sat_xyz[1]], z=[sat_xyz[2]],
        mode="markers", marker=dict(size=7, color="crimson",
                                    line=dict(color="white", width=1)),
        name="satellite",
        hovertemplate=(f"r = {r_sat:.1f} km<br>altitude = {altitude:.1f} km<br>"
                       "x=%{x:.1f}<br>y=%{y:.1f}<br>z=%{z:.1f}<extra>satellite</extra>")))

    # Satellite -> Sun (direction only; the Sun is ~1.5e8 km away)
    L = arrow_len or 1.0 * R
    _arrow(fig, sat_xyz, u * L, "gold", "sat → Sun", width=7)

    # Earth -> Sun direction from the centre, for context
    s_hat = sun_xyz / np.linalg.norm(sun_xyz)
    _arrow(fig, np.zeros(3), s_hat * 2.0 * R, "red", "Earth → Sun (direction)", width=4)

    lim = 2.4 * R
    ax = dict(range=[-lim, lim], showbackground=False, zeroline=False,
              gridcolor="rgba(150,150,150,0.3)")
    fig.update_layout(
        title=(f"altitude {altitude:.1f} km · angle between sat→Sun and "
               f"local zenith: {angle:.1f}°"),
        scene=dict(xaxis=dict(title="x (km)", **ax), yaxis=dict(title="y (km)", **ax),
                   zaxis=dict(title="z (km)", **ax), aspectmode="cube"),
        legend=dict(x=0.01, y=0.99), margin=dict(l=0, r=0, t=50, b=0))
    if show:
        fig.show()
    return fig, u, angle


if __name__ == "__main__":
    sat = [2047.2900751114767, 5124.312255450726, -3978.51714199821]
    sun = 1.496e8 * np.array([0.8, 0.45, 0.4]) / np.linalg.norm([0.8, 0.45, 0.4])
    fig, u, angle = plot_sat_sun(sat, sun, show=False)
    fig.write_html("sat_sun.html")
    print("wrote sat_sun.html", u, angle)