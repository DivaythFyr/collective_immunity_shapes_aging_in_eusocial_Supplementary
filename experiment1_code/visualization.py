import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np
import os

def find_project_root():
    """Find the Bee_eusociality_project folder"""
    current_dir = os.getcwd()
    print(f"Looking for project root starting from: {current_dir}")
    
    # Look for the project folder by going up the directory tree
    path_to_check = current_dir
    
    # Go up the directory tree looking for "Bee_eusociality_project"
    for i in range(10):  # Limit search to 10 levels up
        folder_name = os.path.basename(path_to_check)
        print(f"Checking folder: {folder_name} at path: {path_to_check}")
        
        if folder_name == "Bee_eusociality_project":
            print(f"Found project root: {path_to_check}")
            return path_to_check
        
        parent_path = os.path.dirname(path_to_check)
        
        # If we've reached the root or can't go up further, break
        if parent_path == path_to_check:
            break
            
        path_to_check = parent_path
    
    # If not found, use the assumption that we're in Bee_eusociality_mojo
    fallback_path = os.path.dirname(current_dir)
    print(f"Project root not found, using fallback: {fallback_path}")
    return fallback_path

def load_csv_efficiently(csv_file, max_rows=None, filter_days=None, sample_specific_days=None):
    """
    Load CSV efficiently - either all data for small files, or specific days for large files
    
    Args:
        csv_file: Path to CSV file
        max_rows: Maximum rows to load (for simple row limiting)
        filter_days: Only load data up to this many days
        sample_specific_days: List of specific days to load (for temporal sampling)
    """
    try:
        # Check file size
        file_size_mb = os.path.getsize(csv_file) / (1024 * 1024)
        print(f"✓ File size: {file_size_mb:.1f} MB")
        
        # For small files, load normally
        if file_size_mb < 50 and sample_specific_days is None:
            print("✓ Small file - loading all data...")
            df = pd.read_csv(csv_file)
            if filter_days and 'day' in df.columns:
                df = df[df['day'] <= filter_days].copy()
            return df
        
        # For large files with specific days requested
        if sample_specific_days is not None:
            print(f"✓ Loading specific days: {len(sample_specific_days)} time points")
            print(f"✓ Days to load: {sample_specific_days[:5]}...{sample_specific_days[-5:] if len(sample_specific_days) > 10 else sample_specific_days[5:]}")
            
            # Load data in chunks and filter for specific days
            df_chunks = []
            chunk_size = 50000
            
            for chunk in pd.read_csv(csv_file, chunksize=chunk_size):
                chunk_filtered = chunk[chunk['day'].isin(sample_specific_days)]
                if len(chunk_filtered) > 0:
                    df_chunks.append(chunk_filtered)
            
            if df_chunks:
                df = pd.concat(df_chunks, ignore_index=True)
                print(f"✓ Loaded {df.shape[0]} rows for {len(df['day'].unique())} unique days")
                return df
            else:
                print("✗ No data found for specified days")
                return None
        
        # For large files without specific day requirements - use row limiting
        print(f"⚠️  Large file detected. Using row limiting...")
        limit_rows = min(max_rows or 500000, 500000)
        print(f"✓ Loading first {limit_rows:,} rows...")
        
        df = pd.read_csv(csv_file, nrows=limit_rows)
        if filter_days and 'day' in df.columns:
            df = df[df['day'] <= filter_days].copy()
            
        return df
        
    except Exception as e:
        print(f"✗ Error loading CSV: {e}")
        import traceback
        traceback.print_exc()
        return None

def get_24_temporal_points(min_day, max_day, num_points=24):
    """Calculate exactly 24 evenly distributed time points"""
    if max_day - min_day + 1 <= num_points:
        # If we have fewer days than requested points, return all days
        return list(range(min_day, max_day + 1))
    
    # Calculate exact spacing for num_points
    step_size = (max_day - min_day) / (num_points - 1)
    time_points = []
    
    for i in range(num_points):
        if i == 0:
            day = min_day
        elif i == num_points - 1:
            day = max_day
        else:
            day = min_day + int(round(i * step_size))
        time_points.append(day)
    
    # Remove duplicates while preserving order
    return sorted(list(set(time_points)))

def create_allele_fitness_ranking_plots(allele_csv_file, output_dir, run_id):
    """Create allele fitness ranking line plots using ALL swarms from the final day"""
    try:
        print(f"Creating allele fitness ranking plots...")
        
        # Set matplotlib backend to non-interactive
        plt.switch_backend('Agg')
        
        if not os.path.exists(allele_csv_file):
            print(f"Error: Allele file {allele_csv_file} does not exist!")
            return False
        
        # First, find the final day without loading all data
        print("✓ Finding final day...")
        day_info = pd.read_csv(allele_csv_file, usecols=['day'])
        final_day = day_info['day'].max()
        print(f"✓ Final day (steady state): {final_day}")
        
        # Load ALL data from ONLY the final day
        print(f"✓ Loading ALL swarms from final day {final_day}...")
        df = load_csv_efficiently(allele_csv_file, sample_specific_days=[final_day])
        
        if df is None or df.empty:
            print("✗ No data loaded for final day!")
            return False
        
        # Verify we have final day data
        final_day_data = df[df['day'] == final_day].copy()
        if final_day_data.empty:
            print(f"✗ No data found for final day {final_day}!")
            return False
        
        print(f"✓ Loaded {len(final_day_data)} swarms from final day {final_day}")
        
        # Check data integrity
        required_columns = ['day', 'swarm_id', 'queen_age_allele', 'transition_age_allele', 'forager_age_allele', 'hypersensitivity_allele', 'fitness']
        missing_columns = [col for col in required_columns if col not in final_day_data.columns]
        if missing_columns:
            print(f"Error: Missing columns: {missing_columns}")
            return False
        
        # Sort by fitness (worst to best) and create ranking
        final_day_data = final_day_data.sort_values('fitness', ascending=True).reset_index(drop=True)
        final_day_data['fitness_rank'] = range(len(final_day_data), 0, -1)  # Worst to best rank
        
        print(f"✓ Creating fitness ranking plots for ALL {len(final_day_data)} swarms...")
        print(f"✓ Fitness range: {final_day_data['fitness'].min()} to {final_day_data['fitness'].max()}")
        print(f"✓ Queen age allele range: {final_day_data['queen_age_allele'].min()} to {final_day_data['queen_age_allele'].max()}")
        print(f"✓ Transition age allele range: {final_day_data['transition_age_allele'].min()} to {final_day_data['transition_age_allele'].max()}")
        print(f"✓ Forager age allele range: {final_day_data['forager_age_allele'].min()} to {final_day_data['forager_age_allele'].max()}")
        print(f"✓ hypersensitivity allele range: {final_day_data['hypersensitivity_allele'].min()} to {final_day_data['hypersensitivity_allele'].max()}")
        
        # Create line plots: allele values vs fitness rank
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(18, 6))
        
        # Plot 1: Queen Age Allele vs Fitness Rank
        ax1.plot(final_day_data['queen_age_allele'], final_day_data['fitness_rank'], 'o-', alpha=0.7, markersize=2)
        ax1.set_xlabel('Queen Age Allele Value')
        ax1.set_ylabel(f'Fitness Rank ({len(final_day_data)}=Worst, 1=Best)')
        ax1.set_title('Queen Age Allele vs Fitness Rank')
        ax1.grid(True, alpha=0.3)
        ax1.invert_yaxis()  # Best fitness (1) at top
        
        # Plot 2: Transition Age Allele vs Fitness Rank
        ax2.plot(final_day_data['transition_age_allele'], final_day_data['fitness_rank'], 'o-', alpha=0.7, markersize=2, color='orange')
        ax2.set_xlabel('Transition Age Allele Value')
        ax2.set_ylabel(f'Fitness Rank ({len(final_day_data)}=Worst, 1=Best)')
        ax2.set_title('Transition Age Allele vs Fitness Rank')
        ax2.grid(True, alpha=0.3)
        ax2.invert_yaxis()
        
        # Plot 3: Forager Age Allele vs Fitness Rank
        ax3.plot(final_day_data['forager_age_allele'], final_day_data['fitness_rank'], 'o-', alpha=0.7, markersize=2, color='green')
        ax3.set_xlabel('Forager Age Allele Value')
        ax3.set_ylabel(f'Fitness Rank ({len(final_day_data)}=Worst, 1=Best)')
        ax3.set_title('Forager Age Allele vs Fitness Rank')
        ax3.grid(True, alpha=0.3)
        ax3.invert_yaxis()
        
        # Plot 4: Forager Age Allele vs Fitness Rank
        ax4.plot(final_day_data['hypersensitivity_allele'], final_day_data['fitness_rank'], 'o-', alpha=0.7, markersize=2, color='purple')
        ax4.set_xlabel('Hypersensitivity  Allele Value')
        ax4.set_ylabel(f'Fitness Rank ({len(final_day_data)}=Worst, 1=Best)')
        ax4.set_title('Hypersensitivity Allele vs Fitness Rank')
        ax4.grid(True, alpha=0.3)
        ax4.invert_yaxis()
        
        # Add comprehensive title with stats
        fig.suptitle(f'Allele Values vs Fitness Ranking - Day {final_day} Steady State\n'
                    f'(Run: {run_id}) - All {len(final_day_data)} Swarms', 
                    fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        # Save plots
        ranking_filename = f'{run_id}_allele_fitness_ranking.png'
        ranking_path = os.path.join(output_dir, ranking_filename)
        plt.savefig(ranking_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Fitness ranking plots saved: {ranking_path}")
        print(f"✓ Used ALL {len(final_day_data)} swarms from steady state (day {final_day})")
        return True
        
    except Exception as e:
        print(f"✗ Error creating allele fitness ranking plots: {e}")
        import traceback
        traceback.print_exc()
        return False

def create_allele_tufte_boxplots(allele_csv_file, output_dir, run_id):
    """Create allele temporal distribution Tufte boxplots with exactly 24 time points"""
    try:
        print(f"Creating allele Tufte boxplots...")
        
        # Set matplotlib backend to non-interactive
        plt.switch_backend('Agg')
        
        if not os.path.exists(allele_csv_file):
            print(f"Error: Allele file {allele_csv_file} does not exist!")
            return False
        
        # First, get the time range
        print("✓ Reading time range...")
        day_info = pd.read_csv(allele_csv_file, usecols=['day'])
        min_day = day_info['day'].min()
        max_day = day_info['day'].max()
        
        # Calculate exactly 24 time points
        sample_days = get_24_temporal_points(min_day, max_day, 24)
        print(f"✓ Selected 24 time points: {sample_days}")
        
        # Load only the data we need
        df = load_csv_efficiently(allele_csv_file, sample_specific_days=sample_days)
        if df is None:
            return False
        
        # Check data integrity
        required_columns = ['day', 'swarm_id', 'queen_age_allele', 'transition_age_allele', 'forager_age_allele', 'hypersensitivity_allele', 'fitness']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing columns: {missing_columns}")
            return False
        
        # Prepare data for box plots - ensure we have exactly 24 points
        queen_age_data = []
        transition_age_data = []
        forager_age_data = []
        hypersensitivity_data = []
        box_labels = []
        
        for day in sample_days:
            day_data = df[df['day'] == day]
            if len(day_data) > 0:
                queen_age_data.append(day_data['queen_age_allele'].values)
                transition_age_data.append(day_data['transition_age_allele'].values)
                forager_age_data.append(day_data['forager_age_allele'].values)
                hypersensitivity_data.append(day_data['hypersensitivity_allele'].values)
                
                # Create labels based on simulation length
                if max_day <= 100:
                    box_labels.append(f'Day {day}')
                elif max_day <= 1000:
                    box_labels.append(f'D{day}')
                else:
                    year = day / 365.0
                    if year < 1:
                        box_labels.append(f'D{day}')
                    else:
                        box_labels.append(f'Y{year:.1f}')
            else:
                print(f"Warning: No data for day {day}")
        
        print(f"✓ Successfully prepared data for {len(queen_age_data)} time points")
        
        # Calculate final day statistics
        def calculate_final_statistics(data, allele_name):
            """Calculate min, max, median for the final day"""
            if not data or len(data) == 0:
                return "No data", "No data", "No data"
            
            final_day_data = data[-1]  # Last time point
            if len(final_day_data) == 0:
                return "No data", "No data", "No data"
            
            min_val = np.min(final_day_data)
            max_val = np.max(final_day_data)
            median_val = np.median(final_day_data)
            
            return min_val, max_val, median_val
        
        # Get final day statistics for each allele
        queen_min, queen_max, queen_median = calculate_final_statistics(queen_age_data, "Queen")
        transition_min, transition_max, transition_median = calculate_final_statistics(transition_age_data, "Transition")
        forager_min, forager_max, forager_median = calculate_final_statistics(forager_age_data, "Forager")
        hypersensitivity_min, hypersensitivity_max, hypersensitivity_median = calculate_final_statistics(hypersensitivity_data, "Hypersensitivity")
        
        print(f"✓ Final day statistics:")
        print(f"  Queen: Min={queen_min}, Max={queen_max}, Median={queen_median}")
        print(f"  Transition: Min={transition_min}, Max={transition_max}, Median={transition_median}")
        print(f"  Forager: Min={forager_min}, Max={forager_max}, Median={forager_median}")
        print(f"  Hypersensitivity: Min={hypersensitivity_min}, Max={hypersensitivity_max}, Median={hypersensitivity_median}")
        
        # Calculate dynamic Y-axis limits for each allele type
        def calculate_y_limits(data, padding_factor=0.1):
            """Calculate Y-axis limits based on actual data range with padding"""
            if not data:
                return 0, 100  # Default fallback
            
            # Flatten all data points across all time periods
            all_values = np.concatenate(data) if len(data) > 0 else np.array([])
            
            if len(all_values) == 0:
                return 0, 100
            
            min_val = np.min(all_values)
            max_val = np.max(all_values)
            
            # Add padding (10% of range on each side)
            range_padding = (max_val - min_val) * padding_factor
            y_min = max(0, min_val - range_padding)  # Don't go below 0 for ages
            y_max = max_val + range_padding
            
            return y_min, y_max
        
        # Get Y-axis limits for each allele type
        queen_y_min, queen_y_max = calculate_y_limits(queen_age_data)
        transition_y_min, transition_y_max = calculate_y_limits(transition_age_data)
        forager_y_min, forager_y_max = calculate_y_limits(forager_age_data)
        hypersensitivity_y_min, hypersensitivity_y_max = calculate_y_limits(hypersensitivity_data)
        
        print(f"✓ Dynamic Y-axis limits:")
        print(f"  Queen: {queen_y_min:.1f} to {queen_y_max:.1f}")
        print(f"  Transition: {transition_y_min:.1f} to {transition_y_max:.1f}")
        print(f"  Forager: {forager_y_min:.1f} to {forager_y_max:.1f}")
        print(f"  Hypersensitivity: {hypersensitivity_y_min:.1f} to {hypersensitivity_y_max:.1f}")
        
        # Create the figure with custom styling - FIXED: 5 rows instead of 4
        fig = plt.figure(figsize=(20, 16))  # Slightly taller figure
        gs = plt.GridSpec(5, 1, height_ratios=[3, 3, 3, 3, 1])  # 4 boxplots + 1 stats panel
        
        ax1 = plt.subplot(gs[0])  # Queen boxplot
        ax2 = plt.subplot(gs[1])  # Transition boxplot  
        ax3 = plt.subplot(gs[2])  # Forager boxplot
        ax4 = plt.subplot(gs[3])  # Hypersensitivity boxplot
        ax_stats = plt.subplot(gs[4])  # Statistics panel in separate row
        
        # Define Tufte-style colors
        median_color = '#2c3e50'
        quartile_color = '#34495e'
        outlier_color = '#e74c3c'
        
        def create_tufte_boxplot(ax, data, labels, title, ylabel, y_limits):
            """Create a Tufte-style boxplot with dynamic Y-axis limits"""
            if not data:
                return
            
            # Create regular boxplot but hide most elements
            bp = ax.boxplot(data, labels=labels, patch_artist=False, 
                           showfliers=True, showbox=False, showcaps=False,
                           medianprops={'color': median_color, 'linewidth': 2},
                           whiskerprops={'color': quartile_color, 'linewidth': 1, 'linestyle': '-'},
                           flierprops={'marker': 'o', 'markerfacecolor': outlier_color, 
                                     'markeredgecolor': outlier_color, 'markersize': 3, 'alpha': 0.6})
            
            # Calculate positions and statistics for custom elements
            positions = range(1, len(data) + 1)
            
            for i, pos in enumerate(positions):
                if len(data[i]) > 0:
                    # Calculate quartiles
                    q1 = np.percentile(data[i], 25)
                    q3 = np.percentile(data[i], 75)
                    median = np.percentile(data[i], 50)
                    
                    # Draw quartile lines (horizontal lines at Q1 and Q3)
                    line_width = 0.3
                    ax.hlines(q1, pos - line_width, pos + line_width, colors=quartile_color, linewidth=2)
                    ax.hlines(q3, pos - line_width, pos + line_width, colors=quartile_color, linewidth=2)
                    
                    # Emphasize median
                    ax.hlines(median, pos - line_width, pos + line_width, colors=median_color, linewidth=3)
            
            # Tufte-style formatting with dynamic Y-axis
            ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_ylim(y_limits[0], y_limits[1])  # Dynamic Y-axis limits
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['bottom'].set_linewidth(0.5)
            ax.tick_params(axis='x', rotation=45, labelsize=9)
            ax.tick_params(axis='y', labelsize=10)
            ax.set_facecolor('#fafafa')
        
        # Create Tufte boxplots for each allele type with dynamic Y-axis
        create_tufte_boxplot(ax1, queen_age_data, box_labels, 
                            f'Queen Age Allele Distribution Over Time (Run: {run_id})', 
                            'Queen Age Allele Value', (queen_y_min, queen_y_max))
        
        create_tufte_boxplot(ax2, transition_age_data, box_labels,
                            f'Transition Age Allele Distribution Over Time (Run: {run_id})', 
                            'Transition Age Allele Value', (transition_y_min, transition_y_max))
        
        create_tufte_boxplot(ax3, forager_age_data, box_labels,
                            f'Forager Age Allele Distribution Over Time (Run: {run_id})', 
                            'Forager Age Allele Value', (forager_y_min, forager_y_max))
        
        create_tufte_boxplot(ax4, hypersensitivity_data, box_labels,
                            f'Hypersensitivity Allele Distribution Over Time (Run: {run_id})', 
                            'Hypersensitivity Allele Value', (hypersensitivity_y_min, hypersensitivity_y_max))
        
        # FIXED: Use ax4 instead of undefined 'ax'
        ax4.set_xlabel('Time', fontsize=12, fontweight='bold')
        
        # === CREATE STATISTICS PANEL ===
        ax_stats.axis('off')  # Turn off axes for text panel
        
        # Create statistics text - FIXED: Added missing newline
        stats_text = (
            f'FINAL DAY ALLELE STATISTICS (Day {max_day}):\n\n'
            f'QUEEN AGE:\n'
            f'  Minimum: {queen_min}\n'
            f'  Maximum: {queen_max}\n'
            f'  Median:  {queen_median:.1f}\n\n'
            f'TRANSITION AGE:\n'
            f'  Minimum: {transition_min}\n'
            f'  Maximum: {transition_max}\n'
            f'  Median:  {transition_median:.1f}\n\n'
            f'FORAGER AGE:\n'
            f'  Minimum: {forager_min}\n'
            f'  Maximum: {forager_max}\n'
            f'  Median:  {forager_median:.1f}\n\n'  # Added missing \n here
            f'HYPERSENSITIVITY:\n'
            f'  Minimum: {hypersensitivity_min}\n'
            f'  Maximum: {hypersensitivity_max}\n'
            f'  Median:  {hypersensitivity_median:.1f}'
        )
        
        # Add statistics text to panel
        ax_stats.text(0.02, 0.95, stats_text, transform=ax_stats.transAxes, fontsize=11,
                     verticalalignment='top', horizontalalignment='left',
                     bbox=dict(boxstyle='round', facecolor='#f8f9fa', alpha=0.9,
                             edgecolor='#dee2e6', linewidth=1),
                     fontfamily='monospace')
        
        fig.suptitle('Allele Distribution Over Time (Tufte Style)', fontsize=16, fontweight='bold', y=0.98)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        
        # Save Tufte boxplots
        tufte_filename = f'{run_id}_allele_tufte_boxplots.png'
        tufte_path = os.path.join(output_dir, tufte_filename)
        plt.savefig(tufte_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✓ Tufte boxplots saved: {tufte_path}")
        return True
        
    except Exception as e:
        print(f"✗ Error creating Tufte boxplots: {e}")
        import traceback
        traceback.print_exc()  # Add this for debugging
        return False

def create_fitness_tufte_boxplots(allele_csv_file, output_dir, run_id):
    """Create Tufte-style boxplots for fitness W over time with exactly 24 time points"""
    try:
        print(f"Creating fitness Tufte boxplots...")
        
        # Set matplotlib backend to non-interactive
        plt.switch_backend('Agg')
        
        if not os.path.exists(allele_csv_file):
            print(f"Error: Allele file {allele_csv_file} does not exist!")
            return False
        
        # First, get the time range
        print("✓ Reading time range...")
        day_info = pd.read_csv(allele_csv_file, usecols=['day'])
        min_day = day_info['day'].min()
        max_day = day_info['day'].max()
        
        # Calculate exactly 24 time points
        sample_days = get_24_temporal_points(min_day, max_day, 24)
        print(f"✓ Selected 24 time points: {sample_days}")
        
        # Load only the data we need
        df = load_csv_efficiently(allele_csv_file, sample_specific_days=sample_days)
        if df is None:
            return False
        
        # Check data integrity
        required_columns = ['day', 'swarm_id', 'fitness']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing columns: {missing_columns}")
            return False
        
        # Prepare data for box plots
        fitness_data = []
        box_labels = []
        
        for day in sample_days:
            day_data = df[df['day'] == day]
            if len(day_data) > 0:
                fitness_data.append(day_data['fitness'].values)
                
                # Create labels based on simulation length
                if max_day <= 100:
                    box_labels.append(f'Day {day}')
                elif max_day <= 1000:
                    box_labels.append(f'D{day}')
                else:
                    year = day / 365.0
                    if year < 1:
                        box_labels.append(f'D{day}')
                    else:
                        box_labels.append(f'Y{year:.1f}')
            else:
                print(f"Warning: No data for day {day}")
        
        print(f"✓ Successfully prepared fitness data for {len(fitness_data)} time points")
        
        # Create the figure with custom styling
        fig, ax = plt.subplots(1, 1, figsize=(20, 8))
        
        # Define Tufte-style colors
        median_color = '#2c3e50'
        quartile_color = '#34495e'
        outlier_color = '#e74c3c'
        
        def create_fitness_tufte_boxplot(ax, data, labels, title, ylabel):
            """Create a Tufte-style boxplot for fitness"""
            if not data:
                return
            
            # Filter out empty datasets
            filtered_data = []
            filtered_labels = []
            for i, dataset in enumerate(data):
                if len(dataset) > 0:
                    filtered_data.append(dataset)
                    filtered_labels.append(labels[i])
            
            if len(filtered_data) == 0:
                print(f"Warning: All datasets empty for {title}")
                return
            
            # Create regular boxplot but hide most elements
            bp = ax.boxplot(filtered_data, labels=filtered_labels, patch_artist=False, 
                           showfliers=True, showbox=False, showcaps=False,
                           medianprops={'color': median_color, 'linewidth': 2},
                           whiskerprops={'color': quartile_color, 'linewidth': 1, 'linestyle': '-'},
                           flierprops={'marker': 'o', 'markerfacecolor': outlier_color, 
                                     'markeredgecolor': outlier_color, 'markersize': 3, 'alpha': 0.6})
            
            # Calculate positions and statistics for custom elements
            positions = range(1, len(filtered_data) + 1)
            
            for i, pos in enumerate(positions):
                if len(filtered_data[i]) > 0:
                    # Calculate quartiles
                    q1 = np.percentile(filtered_data[i], 25)
                    q3 = np.percentile(filtered_data[i], 75)
                    median = np.percentile(filtered_data[i], 50)
                    
                    # Draw quartile lines
                    line_width = 0.3
                    ax.hlines(q1, pos - line_width, pos + line_width, colors=quartile_color, linewidth=2)
                    ax.hlines(q3, pos - line_width, pos + line_width, colors=quartile_color, linewidth=2)
                    
                    # Emphasize median
                    ax.hlines(median, pos - line_width, pos + line_width, colors=median_color, linewidth=3)
            
            # Tufte-style formatting
            ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
            ax.set_ylabel(ylabel, fontsize=12)
            ax.set_xlabel('Time', fontsize=12, fontweight='bold')
            ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_linewidth(0.5)
            ax.spines['bottom'].set_linewidth(0.5)
            ax.tick_params(axis='x', rotation=45, labelsize=9)
            ax.tick_params(axis='y', labelsize=10)
            ax.set_facecolor('#fafafa')
        
        # Create fitness Tufte boxplot
        create_fitness_tufte_boxplot(ax, fitness_data, box_labels, 
                                    f'Fitness W Distribution Over Time (Run: {run_id})', 
                                    'Fitness W Value')
        
        plt.tight_layout()
        
        # Save fitness Tufte boxplots
        fitness_filename = f'{run_id}_fitness_tufte_boxplots.png'
        fitness_path = os.path.join(output_dir, fitness_filename)
        plt.savefig(fitness_path, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print(f"✓ Fitness Tufte boxplots saved: {fitness_path}")
        return True
        
    except Exception as e:
        print(f"✗ Error creating fitness Tufte boxplots: {e}")
        import traceback
        traceback.print_exc()
        return False
    
def create_total_population_infection_areaplot(infection_csv_file, output_dir, run_id):
    """Create dual area plots - main plot for large numbers, inset for infected queens"""
    try:
        print(f"Creating total population vs infection by caste area plot...")
        
        # Set matplotlib backend to non-interactive
        plt.switch_backend('Agg')
        
        if not os.path.exists(infection_csv_file):
            print(f"Error: Infection file {infection_csv_file} does not exist!")
            return False
        
        # First, get the time range
        print("✓ Reading time range...")
        day_info = pd.read_csv(infection_csv_file, usecols=['day'])
        min_day = day_info['day'].min()
        max_day = day_info['day'].max()
        
        # Calculate exactly 24 time points
        sample_days = get_24_temporal_points(min_day, max_day, 24)
        print(f"✓ Selected 24 time points: {sample_days}")
        
        # Load both infection and population data
        infection_df = load_csv_efficiently(infection_csv_file, sample_specific_days=sample_days)
        
        # Load population data to get total caste counts
        population_filename = infection_csv_file.replace('_infection_data.csv', '_population_data.csv')
        population_df = load_csv_efficiently(population_filename, sample_specific_days=sample_days)
        
        if infection_df is None or population_df is None:
            return False
        
        # Merge infection and population data
        df = pd.merge(infection_df, population_df, on=['day', 'swarm_id'], how='inner')
        
        # Check data integrity
        required_columns = ['day', 'swarm_id', 'infected_queens', 'infected_nurses', 'infected_foragers', 'queens', 'nurses', 'foragers']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing columns: {missing_columns}")
            return False
        
        # Calculate totals for each day across all swarms
        daily_totals = []
        
        for day in sample_days:
            day_data = df[df['day'] == day]
            if len(day_data) > 0:
                # Total population across all swarms for this day
                total_queens = day_data['queens'].sum()
                total_nurses = day_data['nurses'].sum()
                total_foragers = day_data['foragers'].sum()
                total_bees = total_queens + total_nurses + total_foragers
                
                # Total infected by caste across all swarms for this day
                total_infected_queens = day_data['infected_queens'].sum()
                total_infected_nurses = day_data['infected_nurses'].sum()
                total_infected_foragers = day_data['infected_foragers'].sum()
                
                # Calculate healthy bees by caste
                total_healthy_nurses = total_nurses - total_infected_nurses
                total_healthy_foragers = total_foragers - total_infected_foragers
                total_healthy_queens = total_queens - total_infected_queens
                
                daily_totals.append({
                    'day': day,
                    'total_bees': total_bees,
                    'infected_queens': total_infected_queens,
                    'infected_nurses': total_infected_nurses,
                    'infected_foragers': total_infected_foragers,
                    'healthy_nurses': total_healthy_nurses,
                    'healthy_foragers': total_healthy_foragers,
                    'healthy_queens': total_healthy_queens
                })
            else:
                print(f"Warning: No data for day {day}")
        
        if not daily_totals:
            print("Error: No valid data found for any day")
            return False
        
        # Convert to DataFrame for easier plotting
        plot_df = pd.DataFrame(daily_totals)
        
        print(f"✓ Successfully prepared data for {len(plot_df)} time points")
        print(f"✓ Total bee range: {plot_df['total_bees'].min():,} to {plot_df['total_bees'].max():,}")
        print(f"✓ Healthy nurses range: {plot_df['healthy_nurses'].min():,} to {plot_df['healthy_nurses'].max():,}")
        print(f"✓ Healthy foragers range: {plot_df['healthy_foragers'].min():,} to {plot_df['healthy_foragers'].max():,}")
        print(f"✓ Infected nurses range: {plot_df['infected_nurses'].min():,} to {plot_df['infected_nurses'].max():,}")
        print(f"✓ Infected foragers range: {plot_df['infected_foragers'].min():,} to {plot_df['infected_foragers'].max():,}")
        
        # Create time labels based on simulation length
        if max_day <= 100:
            time_labels = [f'Day {day}' for day in plot_df['day']]
        elif max_day <= 1000:
            time_labels = [f'D{day}' for day in plot_df['day']]
        else:
            time_labels = []
            for day in plot_df['day']:
                year = day / 365.0
                if year < 1:
                    time_labels.append(f'D{day}')
                else:
                    time_labels.append(f'Y{year:.1f}')
        
        # === CREATE DUAL PLOT LAYOUT ===
        fig, (ax_main, ax_inset) = plt.subplots(2, 1, figsize=(16, 12), height_ratios=[3, 1])
        
        x_positions = range(len(plot_df))
        
        # === MAIN PLOT: Detailed breakdown of healthy and infected bees ===
        # Build stacked areas from bottom to top
        
        # Layer 1: Healthy foragers (bottom layer)
        ax_main.fill_between(x_positions, 0, plot_df['healthy_foragers'], 
                           color='lightgreen', alpha=0.7, label='Healthy Foragers')
        
        # Layer 2: Healthy nurses (on top of healthy foragers)
        bottom_healthy_nurses = plot_df['healthy_foragers']
        top_healthy_nurses = bottom_healthy_nurses + plot_df['healthy_nurses']
        ax_main.fill_between(x_positions, bottom_healthy_nurses, top_healthy_nurses,
                           color='lightblue', alpha=0.7, label='Healthy Nurses')
        
        # Layer 3: Infected nurses (on top of healthy nurses)
        bottom_infected_nurses = top_healthy_nurses
        top_infected_nurses = bottom_infected_nurses + plot_df['infected_nurses']
        ax_main.fill_between(x_positions, bottom_infected_nurses, top_infected_nurses,
                           color='orange', alpha=0.8, label='Infected Nurses')
        
        # Layer 4: Infected foragers (on top of infected nurses)
        bottom_infected_foragers = top_infected_nurses
        top_infected_foragers = bottom_infected_foragers + plot_df['infected_foragers']
        ax_main.fill_between(x_positions, bottom_infected_foragers, top_infected_foragers,
                           color='red', alpha=0.8, label='Infected Foragers')
        
        # Add total population outline
        ax_main.plot(x_positions, plot_df['total_bees'], color='black', linewidth=2, alpha=0.9, 
                    label='Total Population Outline')
        
        # Format main plot
        ax_main.set_ylabel('Number of Bees', fontsize=14, fontweight='bold')
        ax_main.set_title(f'Bee Population Composition: Healthy vs Infected by Caste\n(Run: {run_id})', 
                         fontsize=16, fontweight='bold', pad=20)
        
        # Set Y-axis limits for main plot
        max_y = 2001 * 300
        ax_main.set_ylim(0, max_y)
        ax_main.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{x:,.0f}'))
        
        # Add legend for main plot
        ax_main.legend(loc='upper right', fontsize=12, framealpha=0.9)
        ax_main.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        # Set x-axis labels for main plot (every 3rd label to avoid crowding)
        step = max(1, len(x_positions) // 8)  # Show about 8 labels
        ax_main.set_xticks([x_positions[i] for i in range(0, len(x_positions), step)])
        ax_main.set_xticklabels([time_labels[i] for i in range(0, len(time_labels), step)], rotation=45)
        
        # === INSET PLOT: Infected queens only ===
        ax_inset.fill_between(x_positions, 0, plot_df['infected_queens'], 
                            color='purple', alpha=0.8, label='Infected Queens')
        
        # Format inset plot
        ax_inset.set_ylabel('Infected\nQueens', fontsize=12, fontweight='bold', color='purple')
        ax_inset.set_xlabel('Time', fontsize=14, fontweight='bold')
        ax_inset.set_title('Infected Queens (Scale: 0-300 max)', fontsize=12, color='purple', pad=10)
        
        # Set appropriate scale for queens (0 to max 300, but use actual max for better visualization)
        queen_max = max(300, plot_df['infected_queens'].max() * 1.1)
        ax_inset.set_ylim(0, queen_max)
        ax_inset.tick_params(axis='y', labelcolor='purple', labelsize=10)
        ax_inset.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        # Set x-axis labels for inset plot
        ax_inset.set_xticks([x_positions[i] for i in range(0, len(x_positions), step)])
        ax_inset.set_xticklabels([time_labels[i] for i in range(0, len(time_labels), step)], rotation=45, ha='right')
        
        # Add small legend for inset
        ax_inset.legend(loc='upper right', fontsize=10, framealpha=0.9)
        
        # === ADD COMPREHENSIVE STATISTICS BOX ===
        final_total = plot_df['total_bees'].iloc[-1]
        final_infected_queens = plot_df['infected_queens'].iloc[-1]
        final_infected_nurses = plot_df['infected_nurses'].iloc[-1]
        final_infected_foragers = plot_df['infected_foragers'].iloc[-1]
        final_healthy_nurses = plot_df['healthy_nurses'].iloc[-1]
        final_healthy_foragers = plot_df['healthy_foragers'].iloc[-1]
        final_total_infected = final_infected_queens + final_infected_nurses + final_infected_foragers
        final_healthy = final_healthy_nurses + final_healthy_foragers + plot_df['healthy_queens'].iloc[-1]
        final_infection_rate = final_total_infected / final_total if final_total > 0 else 0
        
        # Calculate queen infection rate (out of max 300 queens across all swarms)
        queen_infection_rate = final_infected_queens / 300 if 300 > 0 else 0
        
        stats_text = f'Final Composition:\n' \
                    f'Total Bees: {final_total:,}\n' \
                    f'Healthy Nurses: {final_healthy_nurses:,}\n' \
                    f'Healthy Foragers: {final_healthy_foragers:,}\n' \
                    f'Infected Queens: {final_infected_queens:,} ({queen_infection_rate:.1%} of max 300)\n' \
                    f'Infected Nurses: {final_infected_nurses:,}\n' \
                    f'Infected Foragers: {final_infected_foragers:,}\n' \
                    f'Total Infected: {final_total_infected:,}\n' \
                    f'Overall Infection Rate: {final_infection_rate:.1%}'
        
        ax_main.text(0.02, 0.02, stats_text, transform=ax_main.transAxes, fontsize=10,
                    verticalalignment='bottom', horizontalalignment='left',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))
        
        # === TIGHT LAYOUT AND SAVE ===
        plt.tight_layout(pad=2.0)  # Add some padding between subplots
        
        # Save plot
        areaplot_filename = f'{run_id}_total_population_infection_areaplot.png'
        areaplot_path = os.path.join(output_dir, areaplot_filename)
        plt.savefig(areaplot_path, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Dual area plot saved: {areaplot_path}")
        print(f"✓ Main plot shows detailed caste breakdown (scale: 0-{max_y:,})")
        print(f"✓ Inset plot shows infected queens (scale: 0-{queen_max:.0f})")
        return True
        
    except Exception as e:
        print(f"✗ Error creating dual area plot: {e}")
        import traceback
        traceback.print_exc()
        return False    
    