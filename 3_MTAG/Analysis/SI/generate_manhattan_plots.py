#!/usr/bin/env python3
"""
Generate Manhattan plots for loneliness (trait 3) from MTAG results
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from pathlib import Path

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300

def load_mtag_results(file_path):
    """Load MTAG results file"""
    print(f"Loading {file_path}...")
    df = pd.read_csv(file_path, sep='\t')
    return df

def create_manhattan_plot(df, title, output_path, sig_threshold=5e-8):
    """Create Manhattan plot"""
    # Calculate -log10(p-value)
    df['neglog10p'] = -np.log10(df['mtag_pval'])
    
    # Create chromosome positions for plotting
    df['chr'] = df['CHR'].astype(int)
    df['pos'] = df['BP'].astype(int)
    
    # Sort by chromosome and position
    df = df.sort_values(['chr', 'pos'])
    
    # Calculate cumulative position for x-axis
    df['ind'] = range(len(df))
    df_grouped = df.groupby('chr')
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(16, 6))
    
    # Alternate colors for chromosomes
    colors = ['#1f77b4', '#ff7f0e']
    x_labels = []
    x_labels_pos = []
    
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='neglog10p', 
                  color=colors[num % 2], ax=ax, s=5, alpha=0.6)
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)
    
    # Add significance threshold line
    sig_line = -np.log10(sig_threshold)
    ax.axhline(y=sig_line, color='r', linestyle='--', linewidth=1, 
               label=f'p = {sig_threshold}')
    
    # Formatting
    ax.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(df)])
    ax.set_ylim([0, df['neglog10p'].max() * 1.05])
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    # Add statistics to plot
    sig_count = (df['mtag_pval'] < sig_threshold).sum()
    stats_text = f'Significant SNPs (p < {sig_threshold}): {sig_count}\nTotal SNPs: {len(df):,}'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
            fontsize=10, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved Manhattan plot to {output_path}")
    plt.close()
    
    return sig_count

def create_comparison_plot(df_female, df_male, output_path, sig_threshold=5e-8):
    """Create side-by-side Manhattan plots for male vs female comparison"""
    # Calculate -log10(p-value)
    df_female['neglog10p'] = -np.log10(df_female['mtag_pval'])
    df_male['neglog10p'] = -np.log10(df_male['mtag_pval'])
    
    # Ensure data is sorted
    df_female = df_female.sort_values(['CHR', 'BP'])
    df_male = df_male.sort_values(['CHR', 'BP'])
    
    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 10), sharex=False)
    
    # Calculate positions
    for df, ax, title, color_scheme in [
        (df_female, ax1, 'Females (EUR_Female_MM)', ['#d62728', '#e377c2']),
        (df_male, ax2, 'Males (EUR_Male_MM)', ['#2ca02c', '#98df8a'])
    ]:
        df['chr'] = df['CHR'].astype(int)
        df['ind'] = range(len(df))
        df_grouped = df.groupby('chr')
        
        x_labels = []
        x_labels_pos = []
        
        for num, (name, group) in enumerate(df_grouped):
            group.plot(kind='scatter', x='ind', y='neglog10p', 
                      color=color_scheme[num % 2], ax=ax, s=5, alpha=0.6)
            x_labels.append(name)
            x_labels_pos.append((group['ind'].iloc[-1] + group['ind'].iloc[0]) / 2)
        
        # Add significance threshold line
        sig_line = -np.log10(sig_threshold)
        ax.axhline(y=sig_line, color='r', linestyle='--', linewidth=1, 
                   label=f'p = {sig_threshold}')
        
        # Formatting
        ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xticks(x_labels_pos)
        ax.set_xticklabels(x_labels)
        ax.set_xlim([0, len(df)])
        ax.set_ylim([0, max(df['neglog10p'].max(), 10) * 1.05])
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        sig_count = (df['mtag_pval'] < sig_threshold).sum()
        stats_text = f'Significant SNPs: {sig_count}'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, 
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax2.set_xlabel('Chromosome', fontsize=12, fontweight='bold')
    
    fig.suptitle('Loneliness (Trait 3) MTAG Results: Males vs Females Comparison', 
                 fontsize=16, fontweight='bold', y=1.00)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot to {output_path}")
    plt.close()

def main():
    # Base paths
    base_path = Path('/Users/mahmoudabdelmoneum/Desktop/MIT/UROPS/TsaiKellis_Reorganized/UKBB/GWAS/BOLT-LMM_MTAG')
    mtag_path = base_path / 'MTAG' / 'SI'
    output_base = base_path / 'MTAG_Analysis' / 'SI'
    
    # Stratifications
    stratifications = ['EUR_Female_MM', 'EUR_Male_MM', 'EUR_MM']
    
    print("="*80)
    print("Generating Manhattan plots for Loneliness (Trait 3)")
    print("="*80)
    print()
    
    # Store data for comparison plot
    df_female = None
    df_male = None
    
    # Generate individual Manhattan plots
    for strat in stratifications:
        print(f"\nProcessing {strat}...")
        
        # Load data
        input_file = mtag_path / strat / 'results' / f'SI_{strat}_Output_trait_3.txt'
        output_dir = output_base / strat
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f'manhattan_loneliness_{strat}.png'
        
        # Load and plot
        df = load_mtag_results(input_file)
        
        # Store for comparison
        if strat == 'EUR_Female_MM':
            df_female = df.copy()
        elif strat == 'EUR_Male_MM':
            df_male = df.copy()
        
        # Create plot
        title = f'Loneliness (Trait 3) - MTAG Results\n{strat.replace("_", " ")}'
        sig_count = create_manhattan_plot(df, title, output_file)
        
        print(f"  Found {sig_count} significant SNPs")
        print(f"  Saved to: {output_file}")
    
    # Create comparison plot
    print("\n" + "="*80)
    print("Creating Male vs Female Comparison Plot")
    print("="*80)
    
    comparison_output = output_base / 'loneliness_male_vs_female_comparison.png'
    create_comparison_plot(df_female, df_male, comparison_output)
    
    print("\n" + "="*80)
    print("All plots generated successfully!")
    print("="*80)
    print(f"\nPlots saved in: {output_base}")

if __name__ == '__main__':
    main()



