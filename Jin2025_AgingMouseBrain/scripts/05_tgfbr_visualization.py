"""
============================================================================
TGFBR Pathway: Comprehensive Visualization
Tanycyte ↔ Microglia TGF-beta Superfamily Signaling in Aging
============================================================================
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Arc, ConnectionPatch
import matplotlib.patches as mpatches
import os

RESULT_DIR = "d:/vscode/Jin2025_AgingMouseBrain/results/cellchat/tgfbr_pathway"
OUT_DIR = RESULT_DIR
os.makedirs(OUT_DIR, exist_ok=True)

# ===========================================================================
# 1. Load data
# ===========================================================================
print("[1] Loading TGFBR data...")
tgf_diff = pd.read_csv(os.path.join(RESULT_DIR, "tgfbr_Tany_MG_differential.csv"))

# Sort by max probability
tgf_diff = tgf_diff.sort_values('max_prob', ascending=False)
print(f"  {len(tgf_diff)} TGFBR-related LR pairs between Tany and MG")

# ===========================================================================
# 2. Figure 3: Comprehensive TGFBR pathway summary
# ===========================================================================
print("\n[2] Creating comprehensive TGFBR visualization...")

fig = plt.figure(figsize=(24, 16))

# ---- Panel A: Tany->MG LR pair comparison (Young vs Aged) ----
ax1 = fig.add_subplot(2, 3, 1)

# Prepare data - all TGFBR LRs with non-zero prob
plot_data = tgf_diff.copy()
plot_data = plot_data[plot_data['max_prob'] > 0].head(12)

lr_labels = [l.replace('_', '\n') for l in plot_data['LR_pair']]
y_pos = np.arange(len(plot_data))
height = 0.35

bars1 = ax1.barh(y_pos - height/2, plot_data['Tany_to_MG_young'], height,
                 color='#4DBBD5', label='Young(2m)', edgecolor='white', linewidth=0.5)
bars2 = ax1.barh(y_pos + height/2, plot_data['Tany_to_MG_aged'], height,
                 color='#E64B35', label='Aged(18m)', edgecolor='white', linewidth=0.5)

ax1.set_yticks(y_pos)
ax1.set_yticklabels(lr_labels, fontsize=7)
ax1.set_xlabel('Communication Probability', fontsize=10)
ax1.set_title('A. Tanycyte -> Microglia: TGFBR LR Pairs\n(Young vs Aged)', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)
ax1.set_xlim(0, max(plot_data['Tany_to_MG_young'].max(), plot_data['Tany_to_MG_aged'].max()) * 1.3)

# ---- Panel B: Log2FC bar chart ----
ax2 = fig.add_subplot(2, 3, 2)

# Sort by log2FC
fc_data = tgf_diff[tgf_diff['max_prob'] > 0].copy()
fc_data = fc_data.sort_values('Tany_to_MG_log2FC')
fc_data = fc_data[fc_data['Tany_to_MG_log2FC'].abs() > 0.01].head(12)

colors = ['#E64B35' if x > 0 else '#4DBBD5' for x in fc_data['Tany_to_MG_log2FC']]
lr_labels_short = [l[:45] for l in fc_data['LR_pair']]

ax2.barh(np.arange(len(fc_data)), fc_data['Tany_to_MG_log2FC'], color=colors,
         edgecolor='white', linewidth=0.5)
ax2.set_yticks(np.arange(len(fc_data)))
ax2.set_yticklabels(lr_labels_short, fontsize=7)
ax2.axvline(x=0, color='black', linewidth=0.8)
ax2.set_xlabel('log2 Fold Change (Aged/Young)', fontsize=10)
ax2.set_title('B. TGFBR Changes with Aging\n(Red=Up, Blue=Down)', fontsize=12, fontweight='bold')

# ---- Panel C: TGF-beta Ligand Usage Shift ----
ax3 = fig.add_subplot(2, 3, 3)

# Aggregate by ligand
tgf_diff_clean = tgf_diff.copy()
tgf_diff_clean['ligand'] = tgf_diff_clean['LR_pair'].str.extract(r'^(TGFB[123]|GDF11|INHBA|INHBB|BMP[2467]|NODAL)')[0]
tgf_diff_clean = tgf_diff_clean.dropna(subset=['ligand'])

ligand_summary = tgf_diff_clean.groupby('ligand').agg(
    young_total=('Tany_to_MG_young', 'sum'),
    aged_total=('Tany_to_MG_aged', 'sum')
).reset_index()
ligand_summary['change'] = ligand_summary['aged_total'] - ligand_summary['young_total']
ligand_summary = ligand_summary.sort_values('aged_total', ascending=True)

x_pos = np.arange(len(ligand_summary))
h = 0.35
ax3.barh(x_pos - h/2, ligand_summary['young_total'], h, color='#4DBBD5', label='Young(2m)')
ax3.barh(x_pos + h/2, ligand_summary['aged_total'], h, color='#E64B35', label='Aged(18m)')
ax3.set_yticks(x_pos)
ax3.set_yticklabels(ligand_summary['ligand'], fontsize=9)
ax3.set_xlabel('Total Communication Probability', fontsize=10)
ax3.set_title('C. TGF-beta Ligand Usage: Tany->MG', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)

# ---- Panel D: Global TGFBR Redistribution ----
ax4 = fig.add_subplot(2, 3, 4)

# Load the matrix data
mat_y = pd.read_csv(os.path.join(RESULT_DIR, "tgfbr_matrix_young.csv"), index_col=0)
mat_a = pd.read_csv(os.path.join(RESULT_DIR, "tgfbr_matrix_aged.csv"), index_col=0)
mat_diff = mat_a - mat_y

# MG incoming TGFBR from all sources
mg_incoming_y = mat_y.loc[:, 'Microglia'].drop('Microglia')
mg_incoming_a = mat_a.loc[:, 'Microglia'].drop('Microglia')
mg_incoming_diff = mg_incoming_a - mg_incoming_y
mg_incoming_diff = mg_incoming_diff.sort_values(ascending=True)

# Top sources
top_sources = mg_incoming_diff[mg_incoming_diff.abs() > 1e-12].tail(10)
if len(top_sources) < 5:
    top_sources = mg_incoming_diff.tail(10)

x_pos = np.arange(len(top_sources))
h = 0.35

# Get young and aged values for these sources
src_young = [mg_incoming_y.get(s, 0) for s in top_sources.index]
src_aged = [mg_incoming_a.get(s, 0) for s in top_sources.index]

ax4.barh(x_pos - h/2, src_young, h, color='#4DBBD5', label='Young(2m)')
ax4.barh(x_pos + h/2, src_aged, h, color='#E64B35', label='Aged(18m)')
ax4.set_yticks(x_pos)
ax4.set_yticklabels(top_sources.index, fontsize=8)
ax4.set_xlabel('Total TGFBR Signaling to MG', fontsize=10)
ax4.set_title('D. Microglia TGFBR Incoming: All Sources\n(Sorted by Change)', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)

# ---- Panel E: Ligand-Receptor-Receptor ternary diagram ----
ax5 = fig.add_subplot(2, 3, 5)

# Create a schematic of TGFB2/3 -> TGFBR1/ACVR1 -> TGFBR2/ACVR2 signaling
# Use a network-style plot
ligands = ['TGFB2\n(↓31% aged)', 'TGFB3\n(NEW aged)', 'GDF11\n(→stable)', 'INHBA\n(↓aged)']
receptors1 = ['TGFBR1', 'ACVR1', 'ACVR1B']
receptors2 = ['TGFBR2', 'ACVR2A', 'ACVR2B']

# Position nodes
lx = np.array([0.1, 0.1, 0.1, 0.1])
ly = np.array([0.8, 0.55, 0.3, 0.1])
r1x = np.array([0.5, 0.5, 0.5])
r1y = np.array([0.75, 0.45, 0.15])
r2x = np.array([0.9, 0.9, 0.9])
r2y = np.array([0.75, 0.45, 0.15])

# Draw edges for specific interactions
# TGFB2 -> TGFBR1+TGFBR2
ax5.annotate('', xy=(0.85, 0.18), xytext=(0.15, 0.78),
            arrowprops=dict(arrowstyle='->', color='#4DBBD5', lw=2, alpha=0.7))
# TGFB2 -> ACVR1+TGFBR1
ax5.annotate('', xy=(0.55, 0.13), xytext=(0.15, 0.76),
            arrowprops=dict(arrowstyle='->', color='#4DBBD5', lw=1.5, alpha=0.5))
# TGFB3 -> TGFBR1+TGFBR2 (new in aged)
ax5.annotate('', xy=(0.85, 0.2), xytext=(0.15, 0.53),
            arrowprops=dict(arrowstyle='->', color='#E64B35', lw=2.5, alpha=0.8,
                          connectionstyle='arc3,rad=0.15'))
# TGFB3 -> ACVR1+TGFBR1
ax5.annotate('', xy=(0.55, 0.15), xytext=(0.15, 0.51),
            arrowprops=dict(arrowstyle='->', color='#E64B35', lw=2, alpha=0.6,
                          connectionstyle='arc3,rad=-0.15'))
# GDF11 -> TGFBR1+ACVR2A
ax5.annotate('', xy=(0.88, 0.72), xytext=(0.15, 0.28),
            arrowprops=dict(arrowstyle='->', color='grey', lw=1.5, alpha=0.5))

# Plot nodes
ax5.scatter(lx, ly, s=300, c=['#4DBBD5', '#E64B35', 'grey', 'lightgrey'],
            edgecolors='black', linewidth=1, zorder=5)
ax5.scatter(r1x, r1y, s=200, c='#FFD700', edgecolors='black', linewidth=1, zorder=5)
ax5.scatter(r2x, r2y, s=200, c='#FF8C00', edgecolors='black', linewidth=1, zorder=5)

for i, label in enumerate(ligands):
    ax5.text(lx[i]-0.02, ly[i], label, ha='right', va='center', fontsize=8, fontweight='bold')
for i, label in enumerate(receptors1):
    ax5.text(r1x[i]+0.01, r1y[i]-0.05, label, ha='left', va='top', fontsize=8, color='#B8860B')
for i, label in enumerate(receptors2):
    ax5.text(r2x[i]+0.01, r2y[i]-0.05, label, ha='left', va='top', fontsize=8, color='#D2691E')

ax5.text(0.5, 0.95, 'Type I Receptors', ha='center', fontsize=9, fontweight='bold', color='#B8860B')
ax5.text(0.9, 0.95, 'Type II Receptors', ha='center', fontsize=9, fontweight='bold', color='#D2691E')
ax5.text(0.1, 0.95, 'Ligands', ha='center', fontsize=9, fontweight='bold')

legend_elements = [
    mpatches.Patch(color='#4DBBD5', label='Down in Aged (TGFB2)'),
    mpatches.Patch(color='#E64B35', label='NEW in Aged (TGFB3)'),
    mpatches.Patch(color='grey', label='Stable (GDF11)'),
]
ax5.legend(handles=legend_elements, loc='lower center', fontsize=8, ncol=3)
ax5.set_xlim(0, 1.05)
ax5.set_ylim(-0.05, 1.05)
ax5.axis('off')
ax5.set_title('E. TGF-beta Ligand-Receptor Wiring Diagram\n(Tany -> MG)', fontsize=12, fontweight='bold')

# ---- Panel F: Key Biology Summary ----
ax6 = fig.add_subplot(2, 3, 6)
ax6.axis('off')

summary_text = """F. TGFBR Pathway: Key Biological Insights

1. ALL TGFBR SIGNALING IS TANYCYTE -> MICROGLIA (ONE-WAY)
   Microglia do NOT send TGF-beta signals back to Tanycytes.
   This is a directional homeostatic regulation pathway.

2. TGFB2 IS THE DOMINANT LIGAND BUT DECLINES WITH AGING
   - TGFB2-TGFBR1-TGFBR2: -31% (log2FC=-0.54)
   - TGFB2-ACVR1-TGFBR1: -35% (log2FC=-0.62)
   Loss of TGFB2 signaling => reduced anti-inflammatory tone

3. TGFB3 APPEARS DE NOVO IN AGING (COMPENSATORY?)
   - TGFB3-TGFBR1-TGFBR2: log2FC=+67.7 (absent in young)
   - TGFB3-ACVR1-TGFBR1: log2FC=+66.5 (absent in young)
   Potential compensatory mechanism for TGFB2 loss

4. GDF11 (REJUVENATION FACTOR) IS STABLE
   GDF11-TGFBR1-ACVR2A/B signaling preserved in aging
   Contrasts with TGFB2 decline

5. GLOBAL TGFBR REDISTRIBUTION TO MICROGLIA
   Microglia becomes #1 TGFBR receiver in aged brain
   MG autocrine TGFBR: MASSIVE increase (+2.3e-08)
   Endothelial->MG: +7.2e-09
   SMC, OPC->MG: +2.5e-09 each
   Net effect: MG TGFBR tone maintained but source shifted

6. BIOLOGICAL IMPLICATION
   Tanycyte-derived TGFB2 anti-inflammatory signaling to
   microglia is compromised in aging. This is partially
   compensated by: (a) TGFB3 emergence, (b) increased
   autocrine MG-MG signaling, (c) vascular niche
   (Endothelial/SMC) contribution. However, autocrine TGFB
   in microglia may have different downstream effects than
   paracrine Tanycyte-derived TGFB signaling."""

ax6.text(0.02, 0.98, summary_text, transform=ax6.transAxes,
         fontsize=8.5, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

plt.suptitle('TGFBR Pathway: Tanycyte-Microglia Signaling in Hypothalamic Aging\n'
             'Jin et al. 2025 Nature Dataset - CellChat Analysis',
             fontsize=16, fontweight='bold', y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "03_tgfbr_comprehensive_summary.png"),
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("  [OK] 03_tgfbr_comprehensive_summary.png")

# ===========================================================================
# 3. Figure 4: MG TGFBR source redistribution heatmap
# ===========================================================================
print("\n[3] Creating TGFBR redistribution heatmap...")

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

# Get MG incoming from all sources, top contributors
all_sources = ['Microglia', 'Endothelial', 'SMC', 'OPC', 'VLMC',
               'Glutamatergic_Neuron', 'Pericyte', 'Astrocyte', 'MOL',
               'Tanycyte', 'GABAergic_Neuron', 'Ependymal', 'COP', 'MFOL']

# Use actual data from matrices
available_sources = [s for s in all_sources if s in mat_y.index and s != 'Microglia']
# Add Microglia separately for autocrine
if 'Microglia' in mat_y.index:
    available_sources_full = ['Microglia'] + available_sources
else:
    available_sources_full = available_sources

# Get values
mg_in_y = [mat_y.loc[s, 'Microglia'] if s in mat_y.index else 0 for s in available_sources_full]
mg_in_a = [mat_a.loc[s, 'Microglia'] if s in mat_a.index else 0 for s in available_sources_full]

# Sort by aged values
sort_idx = np.argsort(mg_in_a)
available_sorted = [available_sources_full[i] for i in sort_idx]
mg_in_y_sorted = [mg_in_y[i] for i in sort_idx]
mg_in_a_sorted = [mg_in_a[i] for i in sort_idx]

y = np.arange(len(available_sorted))
h = 0.35

# Panel A: Bar chart
ax = axes[0]
ax.barh(y - h/2, mg_in_y_sorted, h, color='#4DBBD5', label='Young(2m)')
ax.barh(y + h/2, mg_in_a_sorted, h, color='#E64B35', label='Aged(18m)')
ax.set_yticks(y)
ax.set_yticklabels(available_sorted, fontsize=9)
ax.set_xlabel('Total TGFBR Signal to Microglia')
ax.set_title('A. Microglia TGFBR Incoming\n(All Sources)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

# Panel B: Log2 fold change
ax = axes[1]
log2fc_vals = []
for s in available_sorted:
    y_val = mat_y.loc[s, 'Microglia'] if s in mat_y.index else 0
    a_val = mat_a.loc[s, 'Microglia'] if s in mat_a.index else 0
    if y_val > 0:
        log2fc_vals.append(np.log2((a_val + 1e-30) / (y_val + 1e-30)))
    elif a_val > 0:
        log2fc_vals.append(10)  # arbitrary large
    else:
        log2fc_vals.append(0)

colors_fc = ['#E64B35' if x > 0 else '#4DBBD5' for x in log2fc_vals]
ax.barh(y, log2fc_vals, color=colors_fc, edgecolor='white', linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(available_sorted, fontsize=9)
ax.axvline(x=0, color='black', linewidth=0.8)
ax.set_xlabel('log2FC (Aged/Young)')
ax.set_title('B. Change in TGFBR Signal to MG', fontsize=12, fontweight='bold')

# Panel C: Proportion pie-style comparison
ax = axes[2]
# Calculate proportion
young_total = sum(mg_in_y_sorted)
aged_total = sum(mg_in_a_sorted)

young_pct = [100 * v / young_total if young_total > 0 else 0 for v in mg_in_y_sorted]
aged_pct = [100 * v / aged_total if aged_total > 0 else 0 for v in mg_in_a_sorted]

# Focus on top contributors
top_n = 5
# Sort by aged_pct
top_idx = np.argsort(aged_pct)[-top_n:]
others_y = sum(young_pct[i] for i in range(len(young_pct)) if i not in top_idx)
others_a = sum(aged_pct[i] for i in range(len(aged_pct)) if i not in top_idx)

top_sources_names = [available_sorted[i] for i in top_idx]
top_young_pct = [young_pct[i] for i in top_idx]
top_aged_pct = [aged_pct[i] for i in top_idx]

x_pos = np.arange(top_n + 1)
w = 0.35
ax.bar(x_pos[:-1] - w/2, top_young_pct, w, color='#4DBBD5', label='Young(2m)')
ax.bar(x_pos[:-1] + w/2, top_aged_pct, w, color='#E64B35', label='Aged(18m)')
ax.bar(x_pos[-1] - w/2, others_y, w, color='#4DBBD5')
ax.bar(x_pos[-1] + w/2, others_a, w, color='#E64B35')

ax.set_xticks(x_pos)
ax.set_xticklabels(top_sources_names + ['Others'], fontsize=8, rotation=30)
ax.set_ylabel('% of Total TGFBR Signal to MG')
ax.set_title('C. TGFBR Source Composition\n(Top 5 + Others)', fontsize=12, fontweight='bold')
ax.legend(fontsize=9)

plt.suptitle('TGFBR Signaling Redistribution to Microglia in Hypothalamic Aging',
             fontsize=14, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "04_tgfbr_redistribution.png"),
            dpi=200, bbox_inches='tight', facecolor='white')
plt.close()
print("  [OK] 04_tgfb_redistribution.png")

print("\nDone! All figures in:", OUT_DIR)
