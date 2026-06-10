import numpy as np

# Helper function to make zeros visible (assigns small height for visibility)
def ensure_visible_zeros(values, vis_height):
    return [vis_height if (v is not None and not np.isnan(v) and v == 0) else v for v in values]

def draw_limit(ax, limit, max_val, label):
    ax.axhline(y=limit, color='red', linestyle='--', linewidth=2, label=f'{label} ({limit})')
    y_max = max(limit * 1.2, max_val * 1.1) if not np.isnan(max_val) else limit * 1.2
    ax.axhspan(limit, y_max, color='red', alpha=0.1)
    ax.set_ylim(0, y_max)