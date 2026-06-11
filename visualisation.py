import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# Helper function to make zeros visible (assigns small height for visibility)
def ensure_visible_zeros(values, vis_height):
    return [vis_height if (v is not None and not np.isnan(v) and v == 0) else v for v in values]

def draw_limit(ax, limit, max_val, label):
    ax.axhline(y=limit, color='red', linestyle='--', linewidth=2, label=f'{label} ({limit})')
    y_max = max(limit * 1.2, max_val * 1.1) if not np.isnan(max_val) else limit * 1.2
    ax.axhspan(limit, y_max, color='red', alpha=0.1)
    ax.set_ylim(0, y_max)
    
def graph_generator(names: list[str],
                    mol_weights: list[float],
                    hb_donors_local: list[int],
                    hb_acceptors_local: list[int],
                    hb_donors_pubchem: list[int | None],
                    hb_acceptors_pubchem: list[int | None],
                    lipophilicities: list[float],
                    lipinsky_results: list[bool],
                    compare: bool = False):
    # 2. Initialize plots
    fig = plt.figure(figsize=(10, 8))
    gs = fig.add_gridspec(2, 3)
    ax0 = fig.add_subplot(gs[:, 2])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    fig.suptitle("Parameter Analysis - Lipinski's Rule of Five", fontsize=16, fontweight='bold')
    
    x = np.arange(len(names))
    width = 0.35
    
    # --- PLOT 0: lipinsky compatibility ---
    ax0.axis('off')
    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor='green', edgecolor='black', label='Compliant ligand'),
        Rectangle((0, 0), 1, 1, facecolor='red', edgecolor='black', label='Non-compliant ligand'),
        Rectangle((0, 0), 1, 1, facecolor='skyblue', edgecolor='black', label='Local calculation'),
        Rectangle((0, 0), 1, 1, facecolor='orange', edgecolor='black', label='PubChem'),
        Line2D([0], [0], color='red', lw=2, linestyle='--', label='Ro5 limit'),
        Rectangle((0, 0), 1, 1, facecolor='red', alpha=0.1, label='Exceeded limit area'),
        Rectangle((0, 0), 1, 1, facecolor='white', edgecolor='black', hatch='///', label='Non-compliant bar hatch'),
    ]
    ligand_handles = [
        Rectangle((0, 0), 1, 1, facecolor='green' if compliant else 'red', edgecolor='black', label=name)
        for name, compliant in zip(names, lipinsky_results)
    ]
    legend1 = ax0.legend(handles=legend_handles, loc='upper center', frameon=True, facecolor='white', edgecolor='black')
    ax0.add_artist(legend1)
    ax0.legend(handles=ligand_handles, loc='lower center', frameon=True, facecolor='white', edgecolor='black', title='Ligand list')
    

    # --- PLOT 1: Molecular Weight (Limit <= 500 Da) ---
    # Zero is unlikely for MW, but for safety (e.g., 5 Da for visibility)
    mw_vis = ensure_visible_zeros(mol_weights, 5.0) 
    bars = ax1.bar(x, mw_vis, color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
    for bar, res in zip(bars, lipinsky_results):
        if not res:
            bar.set_hatch('///')
    ax1.set_title("Molecular Weight", fontsize=12)
    ax1.set_ylabel("Da")
    max_mw = max(mol_weights) if mol_weights else 500
    draw_limit(ax1, 500, max_mw, "Ro5 Limit")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, rotation=45, ha='right')
    ax1.legend()

    # --- PLOT 2: Lipophilicity (Limit <= 5) ---
    # LogP can be zero, plot thin line 0.05
    lipo_vis = ensure_visible_zeros(lipophilicities, 0.05)
    bars = ax2.bar(x, lipo_vis, color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
    for bar, res in zip(bars, lipinsky_results):
        if not res:
            bar.set_hatch('///')
    ax2.set_title("Lipophilicity (LogP)", fontsize=12)
    ax2.set_ylabel("LogP")
    draw_limit(ax2, 5, 8, "Ro5 Limit")
    ax2.set_ylim(-8, 8)
    ax2.set_xticks(x)
    ax2.set_xticklabels(names, rotation=45, ha='right')
    ax2.legend()

    # --- PLOT 3: HB Donors (Limit <= 5) ---
    # Replace zeros with visual height 0.1
    hbd_loc_vis = ensure_visible_zeros(hb_donors_local, 0.1)
    hbd_pub_vis = ensure_visible_zeros(hb_donors_pubchem, 0.1)
    
    
    if compare:
        bars = ax3.bar(x - width/2, hbd_loc_vis, width, label='Local calculation', color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
        bars = ax3.bar(x + width/2, hbd_pub_vis, width, label='PubChem', color=["orange" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
    else:
        bars = ax3.bar(x, hbd_loc_vis, width, label='Local calculation', color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
    ax3.set_title("Hydrogen Bond Donors", fontsize=12)
    
    max_hbd = max(max(hb_donors_local) if hb_donors_local else 0, max([v for v in hb_donors_pubchem if not np.isnan(v)], default=0))
    draw_limit(ax3, 5, max_hbd, "Ro5 Limit")
    ax3.set_xticks(x)
    ax3.set_xticklabels(names, rotation=45, ha='right' if compare else "center")
    ax3.legend()

    # --- PLOT 4: HB Acceptors (Limit <= 10) ---
    # Replace zeros with visual height 0.2
    hba_loc_vis = ensure_visible_zeros(hb_acceptors_local, 0.2)
    hba_pub_vis = ensure_visible_zeros(hb_acceptors_pubchem, 0.2)
    
    if compare:
        bars = ax4.bar(x - width/2, hba_loc_vis, width, label='Local calculation', color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
        bars = ax4.bar(x + width/2, hba_pub_vis, width, label='PubChem', color=["orange" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
    else:
        bars = ax4.bar(x, hba_loc_vis, width, label='Local calculation', color=["skyblue" if res else "red" for res in lipinsky_results], edgecolor='black', linewidth=0.5)
        for bar, res in zip(bars, lipinsky_results):
            if not res:
                bar.set_hatch('///')
    ax4.set_title("Hydrogen Bond Acceptors", fontsize=12)
    
    max_hba = max(max(hb_acceptors_local) if hb_acceptors_local else 0, max([v for v in hb_acceptors_pubchem if not np.isnan(v)], default=0))
    draw_limit(ax4, 10, max_hba, "Ro5 Limit")
    ax4.set_xticks(x)
    ax4.set_xticklabels(names, rotation=45, ha='right')
    ax4.legend()

    # Finalization and display
    plt.tight_layout(rect=(0, 0, 1, 0.96))
    plt.show()