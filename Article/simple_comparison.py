#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Load data
sage_df = pd.read_csv("results/protein_matrix_from_sage.csv", index_col=0)
openms_df = pd.read_csv("results/protein_matrix_from_openms.csv", index_col=0)

print("=== OpenMS vs Sage Deep Learning Comparison ===")
print(f"📊 Sage baseline:        {len(sage_df):,} proteins × {len(sage_df.columns)} samples")
print(f"🧠 OpenMS+AlphaPeptDeep: {len(openms_df):,} proteins × {len(openms_df.columns)} samples")

# Basic statistics
sage_proteins = set(sage_df.index)
openms_proteins = set(openms_df.index)
common_proteins = sage_proteins & openms_proteins

print(f"\n🔍 Overlap Analysis:")
print(f"  Common proteins:  {len(common_proteins):,}")
print(f"  Sage exclusive:   {len(sage_proteins - openms_proteins):,}")
print(f"  OpenMS exclusive: {len(openms_proteins - sage_proteins):,}")

# Performance metrics
improvement_ratio = len(openms_proteins) / len(sage_proteins)
toyota_target = 10329
vs_toyota = len(openms_proteins) / toyota_target

print(f"\n📈 Performance Summary:")
print(f"  Improvement over Sage: {improvement_ratio:.1f}x")
print(f"  vs Toyota et al. target: {vs_toyota:.1f}x ({len(openms_proteins) - toyota_target:+,})")

# Create comparison figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('OpenMS Deep Learning vs Sage Library-Free Comparison', fontweight='bold')

# Detection comparison
ax1 = axes[0,0]
methods = ['Sage\n(Library-free)', 'OpenMS\n(Deep Learning)']
counts = [len(sage_proteins), len(openms_proteins)]
colors = ['#ff7f0e', '#2ca02c']
bars = ax1.bar(methods, counts, color=colors, alpha=0.7)
ax1.set_ylabel('Proteins Detected')
ax1.set_title('Detection Performance')
for bar, count in zip(bars, counts):
    ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 500,
            f'{count:,}', ha='center', va='bottom', fontweight='bold')

# Improvement ratio
ax2 = axes[0,1]
ax2.bar(['Improvement\nRatio'], [improvement_ratio], color='#d62728', alpha=0.7)
ax2.set_ylabel('Fold Change')
ax2.set_title('Detection Improvement')
ax2.text(0, improvement_ratio + 0.3, f'{improvement_ratio:.1f}x',
         ha='center', va='bottom', fontweight='bold', fontsize=14)

# Toyota comparison
ax3 = axes[1,0]
targets = ['Toyota\nSage', 'Toyota\nTarget', 'Our\nOpenMS']
target_values = [2110, toyota_target, len(openms_proteins)]
target_colors = ['#ff7f0e', '#1f77b4', '#2ca02c']
bars3 = ax3.bar(targets, target_values, color=target_colors, alpha=0.7)
ax3.set_ylabel('Proteins')
ax3.set_title('vs Toyota et al. 2025')
for bar, val in zip(bars3, target_values):
    ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 200,
            f'{val:,}', ha='center', va='bottom', fontweight='bold')

# Overlap pie chart
ax4 = axes[1,1]
overlap_data = [
    len(sage_proteins - openms_proteins),
    len(common_proteins),
    len(openms_proteins - sage_proteins)
]
labels = ['Sage Only', 'Common', 'OpenMS Only']
overlap_colors = ['#ff7f0e', '#17becf', '#2ca02c']
ax4.pie(overlap_data, labels=labels, colors=overlap_colors, autopct='%1.1f%%', startangle=90)
ax4.set_title('Protein Overlap')

plt.tight_layout()
plt.savefig('results/figures/openms_sage_comparison.png', dpi=300, bbox_inches='tight')
plt.show()

print(f"\n🎉 Analysis complete! Figure saved: results/figures/openms_sage_comparison.png")

# Summary for article
print(f"\n📝 Article Summary:")
print(f"{'='*60}")
print(f"OpenMS + AlphaPeptDeep achieved {len(openms_proteins):,} protein detections")
print(f"vs Sage baseline of {len(sage_proteins):,} proteins")
print(f"Improvement: {improvement_ratio:.1f}x increase")
print(f"Toyota et al. target achievement: {vs_toyota:.1f}x")
print(f"Technical advantage: Deep learning spectral prediction + MBR")
print(f"{'='*60}")
