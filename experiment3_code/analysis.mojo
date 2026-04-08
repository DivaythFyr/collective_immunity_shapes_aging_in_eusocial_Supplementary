from types import *
from simulation_core import create_swarms_with_ranges, generate_swarm_with_chromosome
from pathogen_system import count_working_workers, is_queen_alive_and_fertile
from python import Python
from math import *

fn count_time_sec(start: UInt, end: UInt) -> Float64:
    var total_time = end - start
    var seconds_total_time = total_time / 1e+9
    return seconds_total_time

fn calculate_fitness(swarm: Swarm, pathogens: List[Pathogen]) -> Int:
    """Calculate fitness W for a swarm."""
    var working_workers = count_working_workers(swarm, pathogens)
    var working_nurses = working_workers[0]
    var working_foragers = working_workers[1]
    var queen_fertile_multiplier = 1 if is_queen_alive_and_fertile(swarm, pathogens) else 0

    # Apply constraints for working nurses and foragers
    # There can't be more than 160 impactful nurses and 40 impactful foragers for fitness
    working_nurses = min(working_nurses, MAX_NURSE_CONTRIBUTION)
    working_foragers = min(working_foragers, MAX_FORAGER_CONTRIBUTION)
    
    # Calculate base fitness: working workers - newborns that hatched today
    var base_fitness = (working_nurses * working_foragers) - (swarm.newborns_today*1000) - (len(swarm.bees) * 30)
    
    # Apply queen fertility multiplier
    var fitness = base_fitness * queen_fertile_multiplier
    
    return fitness

fn weighted_select_parents(swarms: List[Swarm], pathogens: List[Pathogen]) -> Tuple[Int, Int]:
    """Select two different parent swarms based on fitness weights. Returns (parent1_index, parent2_index)."""
    var living_swarms = List[SwarmFitness]()
    var total_fitness = 0
    
    # Get all living swarms and their fitness
    for i in range(len(swarms)):
        if len(swarms[i].bees) > 0:  # Only living swarms
            var fitness = calculate_fitness(swarms[i], pathogens)
            living_swarms.append(SwarmFitness(i, fitness))
            total_fitness += fitness
    
    # If no living swarms or insufficient living swarms, return -1
    if len(living_swarms) < 2:
        return (-1, -1)
    
    # If total fitness is 0, select randomly
    if total_fitness == 0:
        var parent1_idx = Int(random_si64(0, len(living_swarms)))
        var parent2_idx = Int(random_si64(0, len(living_swarms) - 1))
        if parent2_idx >= parent1_idx:
            parent2_idx += 1
        return (living_swarms[parent1_idx].swarm_index, living_swarms[parent2_idx].swarm_index)
    
    # Weighted selection for first parent
    var random_val1 = random_float64() * Float64(total_fitness)
    var cumulative_fitness = 0.0
    var parent1_idx = -1
    
    for i in range(len(living_swarms)):
        cumulative_fitness += Float64(living_swarms[i].fitness)
        if random_val1 <= cumulative_fitness:
            parent1_idx = i
            break
    
    # Remove selected parent from consideration for second parent
    var remaining_swarms = List[SwarmFitness]()
    var remaining_total_fitness = 0
    
    for i in range(len(living_swarms)):
        if i != parent1_idx:
            remaining_swarms.append(living_swarms[i].copy())
            remaining_total_fitness += living_swarms[i].fitness
    
    # Weighted selection for second parent
    var parent2_idx = -1
    if remaining_total_fitness == 0:
        # Random selection if all remaining have 0 fitness
        var random_idx = Int(random_si64(0, len(remaining_swarms)))
        parent2_idx = remaining_swarms[random_idx].swarm_index
    else:
        var random_val2 = random_float64() * Float64(remaining_total_fitness)
        cumulative_fitness = 0.0
        
        for i in range(len(remaining_swarms)):
            cumulative_fitness += Float64(remaining_swarms[i].fitness)
            if random_val2 <= cumulative_fitness:
                parent2_idx = remaining_swarms[i].swarm_index
                break
    
    return (living_swarms[parent1_idx].swarm_index, parent2_idx)

fn count_infected_bees(swarms: List[Swarm], pathogens: List[Pathogen]) -> List[Int]:
    """Count infected bees by pathogen type."""
    var infected_counts = List[Int]()
    
    # Handle empty pathogen list safely
    if len(pathogens) == 0:
        return infected_counts.copy()  # Return empty list immediately
    
    # Initialize counts for each pathogen
    for _ in range(len(pathogens)):
        infected_counts.append(0)
    
    # Count infections
    for swarm in swarms:
        for bee in swarm.bees:
            # Check if bee has any infection
            if bee.infected != 0:
                # Check each pathogen bit
                for pathogen_idx in range(len(pathogens)):
                    var pathogen_bit = UInt8(1) << pathogen_idx
                    # If bee is infected by this pathogen, increment count
                    if (bee.infected & pathogen_bit) != 0:
                        infected_counts[pathogen_idx] += 1
                # Note: A bee can be counted for multiple pathogens now
    
    return infected_counts.copy()

fn count_infected_castes(swarms: List[Swarm]) -> Tuple[Int, Int, Int]:
    """Count all infected castes throughout all swarms."""
    var infected_queens = 0
    var infected_nurses = 0
    var infected_foragers = 0
    
    # Count infections
    for swarm in swarms:
        for bee in swarm.bees:
            if bee.infected != 0:
                if bee.caste_type == QUEEN:
                    infected_queens += 1
                elif bee.caste_type == NURSE:
                    infected_nurses += 1
                elif bee.caste_type == FORAGER:
                    infected_foragers += 1

    return (infected_queens, infected_nurses, infected_foragers)

fn count_healthy_castes(swarms: List[Swarm]) -> Tuple[Int, Int, Int]:
    """Count all healthy castes throughout all swarms."""
    var healthy_queens = 0
    var healthy_nurses = 0
    var healthy_foragers = 0
    
    # Count infections
    for swarm in swarms:
        for bee in swarm.bees:
            if bee.infected == 0:
                if bee.caste_type == QUEEN:
                    healthy_queens += 1
                elif bee.caste_type == NURSE:
                    healthy_nurses += 1
                elif bee.caste_type == FORAGER:
                    healthy_foragers += 1

    return (healthy_queens, healthy_nurses, healthy_foragers)

fn count_alive_swarms(swarms: List[Swarm]) -> Int:
    """Count swarms that are still alive (have at least one bee)."""
    var alive_count = 0
    for swarm in swarms:
        if len(swarm.bees) > 0:
            alive_count += 1
    return alive_count

fn record_swarm_death(mut simulation_data: SimulationData):
    """Record that a swarm has died (increment cumulative counter)."""
    simulation_data.total_swarms_died += 1

fn print_simulation_stats(day: Int, swarms: List[Swarm], pathogens: List[Pathogen], simulation_data: SimulationData):
    """Print statistics for the current day."""
    var total_bees = 0
    var total_queens = 0
    var total_nurses = 0
    var total_foragers = 0
    var alive_swarms = count_alive_swarms(swarms)
    
    for swarm in swarms:
        total_bees += len(swarm.bees)
        var counts = count_bees_by_caste(swarm)
        total_queens += counts[0]
        total_nurses += counts[1]
        total_foragers += counts[2]
    
    var infected_counts = count_infected_bees(swarms, pathogens)
    
    print("Day", day, ": Total bees =", total_bees, 
          "| Queens =", total_queens, 
          "| Nurses =", total_nurses, 
          "| Foragers =", total_foragers,
          "| Alive swarms =", alive_swarms,
          "| Total swarms died =", simulation_data.total_swarms_died)
    
    # Print infections for each pathogen
    print("  Infections:", end="")
    if len(pathogens) > 0:  # ← ADD THIS CHECK
        for i in range(len(pathogens)):
            print(" " + String(pathogens[i].ID) + " =", infected_counts[i], end="")
            if i < len(pathogens) - 1:
                print(" |", end="")
    else:
        print(" None (no pathogens)", end="")
    print()

fn find_project_root() -> String:
    """Find the Bee_eusociality_project folder."""
    try:
        var os_module = Python.import_module("os")
        
        # Get current working directory
        var current_dir = os_module.getcwd()
        
        # Look for the project folder by going up the directory tree
        var path_to_check = current_dir
        
        # Go up the directory tree looking for "Bee_eusociality_project"
        for _ in range(10):  # Limit search to 10 levels up
            var folder_name = os_module.path.basename(path_to_check)
            
            if String(folder_name) == "Bee_eusociality_project":
                return String(path_to_check)
            
            var parent_path = os_module.path.dirname(path_to_check)
            
            # If we've reached the root or can't go up further, break
            if parent_path == path_to_check:
                break
                
            path_to_check = parent_path
        
        # If not found, use the assumption that we're in Bee_eusociality_mojo
        var fallback_path = os_module.path.dirname(current_dir)
        print("Project root not found, using fallback:", fallback_path)
        return String(fallback_path)
        
    except:
        print("Warning: Could not determine project root, using fallback")
        return "../"

fn get_output_directory() -> String:
    """Get the correct output directory path in the project root."""
    try:
        var os_module = Python.import_module("os")
        
        # Find project root
        var project_root = find_project_root()
        
        # Create output path in the project root
        var output_dir = os_module.path.join(project_root, "output")
        
        return String(output_dir)
        
    except:
        print("Warning: Could not determine output directory, using default")
        return "../output"

fn create_output_directory():
    """Create output directory if it doesn't exist."""
    try:
        var os_module = Python.import_module("os")
        var output_dir = get_output_directory()
        
        # Check if output directory exists, create if not
        if not os_module.path.exists(output_dir):
            _ = os_module.makedirs(output_dir)
             
    except:
        print("Warning: Could not create output directory")

fn export_simulation_data_to_csv(simulation_data: SimulationData):
    """Export simulation data to CSV files in output folder."""
    try:
        # Create output directory first
        create_output_directory()
        
        # Get output directory path
        var output_dir = get_output_directory()
        
        print("Starting CSV export to:", output_dir)
        
        var os_module = Python.import_module("os")
        var builtins = Python.import_module("builtins")
        
        # Export population data with streaming
        print("Exporting population data...")
        var csv_filename = simulation_data.run_id + "_population_data.csv"
        var csv_path = os_module.path.join(output_dir, csv_filename)
        
        var file = builtins.open(csv_path, "w")
        _ = file.write("day,swarm_id,queens,nurses,foragers\n")
        
        # Write data incrementally to avoid memory issues
        for day in range(simulation_data.days):
            var day_data = simulation_data.daily_data[day].copy()
            for swarm_id in range(len(day_data)):
                var swarm_data = day_data[swarm_id].copy()
                var line = String(day) + "," + String(swarm_id) + "," + String(swarm_data.queens) + "," + String(swarm_data.nurses) + "," + String(swarm_data.foragers) + "\n"
                _ = file.write(line)
        
        _ = file.close()
        print("Population data exported to:", csv_path)
        
        # Export allele data with streaming
        print("Exporting allele data...")
        var allele_filename = simulation_data.run_id + "_allele_data.csv"
        var allele_path = os_module.path.join(output_dir, allele_filename)
        
        var allele_file = builtins.open(allele_path, "w")
        _ = allele_file.write("day,swarm_id,queen_age_allele,transition_age_allele,forager_age_allele,hypersensitivity_allele,fitness\n")
        
        # Write allele data incrementally
        for day in range(simulation_data.days):
            var day_allele_data = simulation_data.daily_allele_data[day].copy()
            for swarm_idx in range(len(day_allele_data)):
                var allele_data = day_allele_data[swarm_idx].copy()
                var allele_line = String(day) + "," + String(allele_data.swarm_id) + "," + String(allele_data.queen_age_allele) + "," + String(allele_data.transition_age_allele) + "," + String(allele_data.forager_age_allele) + "," + String(allele_data.hypersensitivity_allele) + "," + String(allele_data.fitness) + "\n"
                _ = allele_file.write(allele_line)
        
        _ = allele_file.close()
        print("Allele data exported to:", allele_path)
        
        # Export infection data with streaming
        print("Exporting infection data...")
        var infection_filename = simulation_data.run_id + "_infection_data.csv"
        var infection_path = os_module.path.join(output_dir, infection_filename)
        
        var infection_file = builtins.open(infection_path, "w")
        _ = infection_file.write("day,swarm_id,infected_queens,infected_nurses,infected_foragers\n")
        
        # Write infection data incrementally
        for day in range(simulation_data.days):
            var day_infection_data = simulation_data.daily_infection_data[day].copy()
            for swarm_idx in range(len(day_infection_data)):
                var infection_data = day_infection_data[swarm_idx].copy()
                var infection_line = String(day) + "," + String(swarm_idx) + "," + String(infection_data.infected_queens) + "," + String(infection_data.infected_nurses) + "," + String(infection_data.infected_foragers) + "\n"
                _ = infection_file.write(infection_line)
        
        _ = infection_file.close()
        print("Infection data exported to:", infection_path)
        print("✓ All CSV exports completed successfully!")
        
    except e:
        print("Failed to export simulation data to CSV. Error:", e)

fn create_animated_visualization(simulation_data: SimulationData):
    """Create animated visualization by calling Python script - MODIFIED TO SAVE ONLY TUFTE BOXPLOTS AND FITNESS RANKING."""
    try:
        # COMMENTED OUT - CSV export to save disk space
        print("Exporting simulation data to CSV...")
        export_simulation_data_to_csv(simulation_data)
        print("✓ CSV export completed successfully!")
        
        print("Creating ONLY allele tufte boxplots and fitness ranking plots...")
        
        # Get output directory path
        var output_dir = get_output_directory()
        
        # Add current directory to Python path
        Python.add_to_path(".")
        
        print("✓ Python path configured")
        
        # Import the custom visualization module
        var visualization = Python.import_module("visualization")
        print("✓ Visualization module imported")
        
        # Create full CSV paths with run ID prefix
        var os_module = Python.import_module("os")
        var population_filename = simulation_data.run_id + "_population_data.csv"
        var allele_filename = simulation_data.run_id + "_allele_data.csv"
        var infection_filename = simulation_data.run_id + "_infection_data.csv"
        var population_path = os_module.path.join(output_dir, population_filename)
        var allele_path = os_module.path.join(output_dir, allele_filename)
        var infection_path = os_module.path.join(output_dir, infection_filename)
        
        print("✓ File paths created")
        print("  Allele path:", allele_path)

        # NEW - Fitness Tufte boxplots
        print("Creating fitness Tufte boxplots...")
        try:
            _ = visualization.create_fitness_tufte_boxplots(allele_path, output_dir, simulation_data.run_id)
            print("✓ Fitness Tufte boxplots completed")
        except e:
            print("✗ Fitness Tufte boxplots failed:", e)

        try:
            _ = visualization.create_total_population_infection_areaplot(infection_path, output_dir, simulation_data.run_id)
            print("✓ Total population vs infection area plot completed")
        except e:
            print("✗ Total population vs infection area plot failed:", e)
        
        # KEEP - Allele fitness ranking plots
        print("Creating allele fitness ranking plots...")
        try:
            _ = visualization.create_allele_fitness_ranking_plots(allele_path, output_dir, simulation_data.run_id)
            print("✓ Allele fitness ranking plots completed")
        except e:
            print("✗ Allele fitness ranking failed:", e)
        
        # KEEP - Tufte-style allele temporal boxplots
        print("Creating Tufte-style allele temporal boxplots...")
        try:
            _ = visualization.create_allele_tufte_boxplots(allele_path, output_dir, simulation_data.run_id)
            print("✓ Tufte-style boxplots completed")
        except e:
            print("✗ Tufte-style boxplots failed:", e)
        
        print("🎉 Selected visualizations completed!")

        print("Cleaning up CSV files to save disk space...")
        try:
            if os_module.path.exists(population_path):
                _ = os_module.remove(population_path)
                print("✓ Deleted population CSV")
            if os_module.path.exists(allele_path):
                _ = os_module.remove(allele_path)
                print("✓ Deleted allele CSV")
            if os_module.path.exists(infection_path):
                _ = os_module.remove(infection_path)
                print("✓ Deleted infection CSV")
        except e:
            print("⚠️  Could not delete some CSV files:", e)
        
    except e:
        print("Failed to create visualization. Error:", e)
        print("Make sure you have the visualization.py file and required packages installed:")
        print("pip install matplotlib pandas pillow")

# Check steady state on alleles
fn check_steady_state(simulation_data: SimulationData, lookback_days: Int = 100, 
                     coefficient_variation_threshold: Float64 = 0.05, 
                     mean_difference_threshold: Float64 = 0.01,  # Now 1% difference
                     min_stable_checks: Int = 10, day_check_window: Int = 500) -> Bool:
    """Check if alleles have reached steady state based on coefficient of variation AND mean stability over time."""
    
    # Need enough data to analyze
    if simulation_data.days < lookback_days + min_stable_checks * day_check_window:
        return False
    
    # Check the last 'min_stable_checks' periods (each day_check_window days apart)
    var stable_periods = 0
    var previous_queen_mean: Float64 = -1.0
    var previous_transition_mean: Float64 = -1.0
    var previous_forager_mean: Float64 = -1.0
    # var previous_hypersensitibity_mean: Float64 = -1.0
    
    for check_idx in range(min_stable_checks):
        # Calculate which day range to check (going backwards)
        var end_day = simulation_data.days - 1 - (check_idx * day_check_window)
        var start_day = max(0, end_day - lookback_days)
        
        if start_day >= end_day:
            continue
        
        # Collect allele values for this period
        var queen_age_allele_values = List[Float64]()
        var transition_age_allele_values = List[Float64]()
        var forager_age_allele_values = List[Float64]()
        # var hypersensitivity_allele_values = List[Float64]()
        
        for day in range(start_day, end_day + 1):
            if day < len(simulation_data.daily_allele_data):
                var day_data = simulation_data.daily_allele_data[day].copy()
                for swarm_data in day_data:
                    queen_age_allele_values.append(Float64(swarm_data.queen_age_allele))
                    transition_age_allele_values.append(Float64(swarm_data.transition_age_allele))
                    forager_age_allele_values.append(Float64(swarm_data.forager_age_allele))
                    # hypersensitivity_allele_values.append(Float64(swarm_data.hypersensitivity_allele))
        
        # Calculate coefficient of variation for each allele type
        var queen_coefficient_variation = calculate_coefficient_of_variation(queen_age_allele_values)
        var transition_coefficient_variation = calculate_coefficient_of_variation(transition_age_allele_values)
        var forager_coefficient_variation = calculate_coefficient_of_variation(forager_age_allele_values)
        # var hypersensitivity_coefficient_variation = calculate_coefficient_of_variation(hypersensitivity_allele_values)

        # Calculate means for this period
        var queen_mean = calculate_mean(queen_age_allele_values)
        var transition_mean = calculate_mean(transition_age_allele_values)
        var forager_mean = calculate_mean(forager_age_allele_values)
        # var hypersensitivity_mean = calculate_mean(hypersensitivity_allele_values)
        
        # Check BOTH conditions: low variation AND stable means
        var low_variation = (
                            queen_coefficient_variation <= coefficient_variation_threshold
                            and
                            transition_coefficient_variation <= coefficient_variation_threshold
                            and 
                            forager_coefficient_variation <= coefficient_variation_threshold
                            # and
                            # hypersensitivity_coefficient_variation <= coefficient_variation_threshold
                           )
        
        var stable_means = True
        
        # For the first period, we just record the means (nothing to compare against)
        if check_idx > 0:
            # Calculate percentage differences
            var queen_mean_diff_pct = calculate_percentage_difference(previous_queen_mean, queen_mean)
            var transition_mean_diff_pct = calculate_percentage_difference(previous_transition_mean, transition_mean)
            var forager_mean_diff_pct = calculate_percentage_difference(previous_forager_mean, forager_mean)
            # var hypersensitivity_mean_diff_pct = calculate_percentage_difference(previous_hypersensitibity_mean, hypersensitivity_mean)

            stable_means = (
                        queen_mean_diff_pct <= mean_difference_threshold
                        and
                        transition_mean_diff_pct <= mean_difference_threshold 
                        and 
                        forager_mean_diff_pct <= mean_difference_threshold
                        # and
                        # hypersensitivity_mean_diff_pct <= mean_difference_threshold
                        )
        
        # Update previous means for next comparison
        previous_queen_mean = queen_mean
        previous_transition_mean = transition_mean
        previous_forager_mean = forager_mean
        # previous_hypersensitibity_mean = hypersensitivity_mean
        
        # Check if both conditions are satisfied
        if low_variation and stable_means:
            stable_periods += 1
            print("    ✓ Period is stable (low variation + stable means)")
        else:
            var reasons = List[String]()
            if not low_variation:
                reasons.append("high variation")
            if not stable_means and check_idx > 0:  # Only show mean reason if we had comparison
                reasons.append("unstable means")
            break
    
    var is_steady = stable_periods >= min_stable_checks
    if is_steady:
        print("🎯 STEADY STATE REACHED! All alleles stable for", min_stable_checks, "consecutive checks.")
    else:
        pass
    
    return is_steady

fn calculate_coefficient_of_variation(values: List[Float64]) -> Float64:
    """Calculate coefficient of variation (CV = std_dev / mean). Returns relative variability (0.0 to ~1.0+)."""
    if len(values) == 0:
        return 0.0
    
    # Calculate mean
    var sum_val = 0.0
    for value in values:
        sum_val += value
    var mean = sum_val / Float64(len(values))
    
    # Avoid division by zero
    if mean == 0.0:
        return 0.0
    
    # Calculate variance
    var variance_sum = 0.0
    for value in values:
        var diff = value - mean
        variance_sum += diff * diff
    var variance = variance_sum / Float64(len(values))
    
    # Calculate standard deviation
    var std_dev = sqrt(variance)
    
    # Return coefficient of variation (relative measure)
    return std_dev / mean


fn calculate_mean(values: List[Float64]) -> Float64:
    """Calculate mean, used in steady state."""
    if len(values) == 0:
        return 0.0

    Sum = 0.0
    for value in values:
        Sum = Sum + value
    
    return Sum / len(values)

fn calculate_percentage_difference(x1: Float64, x2: Float64) -> Float64:
    """Calculate percentage difference between two numbers. Especially used between two means in check steady state function."""
    if x1 == 0.0:
        return 0.0  # Avoid division by zero
    
    return abs(x2 - x1) / abs(x1)
    
fn bubble_sort(mut values: List[Float64]) -> List[Float64]:
    """Bubble sort implementation - simple but O(n²)."""
    var n = len(values)
    
    for i in range(n):
        var swapped = False
        for j in range(0, n - i - 1):
            if values[j] > values[j + 1]:
                # Swap elements
                var temp = values[j]
                values[j] = values[j + 1]
                values[j + 1] = temp
                swapped = True
        
        # Early exit if already sorted
        if not swapped:
            break
    
    return values.copy()

fn calculate_median(mut values: List[Float64]) -> Float64:
    """Calculate median after bubble sorting. If uneven number of elements - just the middle value. If even - mean of two central values."""
    var index = len(values) // 2

    values = values.copy()
    var sorted_values = bubble_sort(values)

    var median: Float64
    if len(sorted_values) % 2 ==1:
        median = sorted_values[index]
    else:
        median = (sorted_values[index] + sorted_values[index-1]) / 2

    return median

fn calculate_final_infection_proportion(simulation_data: SimulationData) -> Tuple[Float64, Int, Int]:
    """Calculate infected bee proportion on the last day.
    Returns: (infection_proportion, total_infected, total_bees).
    """
    if simulation_data.days == 0:
        return (0.0, 0, 0)
    
    # Get the last day's infection data
    var last_day_idx = simulation_data.days - 1
    var last_day_infection_data = simulation_data.daily_infection_data[last_day_idx].copy()
    
    # Get the last day's population data
    var last_day_population_data = simulation_data.daily_data[last_day_idx].copy()
    
    var total_infected = 0
    var total_bees = 0
    
    # Sum across all swarms for the last day
    for i in range(len(last_day_infection_data)):
        var infection_data = last_day_infection_data[i].copy()
        var population_data = last_day_population_data[i].copy()
        
        # Count infected bees
        total_infected += (infection_data.infected_queens + 
                          infection_data.infected_nurses + 
                          infection_data.infected_foragers)
        
        # Count total bees
        total_bees += (population_data.queens + 
                      population_data.nurses + 
                      population_data.foragers)
    
    # Calculate proportion
    var infection_proportion = 0.0
    if total_bees > 0:
        infection_proportion = Float64(total_infected) / Float64(total_bees)
    
    return (infection_proportion, total_infected, total_bees)

### Collecting basic day statistics for each day of a run ###
fn collect_day_statistics(simulation_data: SimulationData, swarms: List[Swarm], pathogens: List[Pathogen]):
    """Collect each day statistics in the following order and append to CSV file:.

    day
    fitness_mean
    fitness_median
    fitness_CV
    hypersensitivity_mean
    hypersensitivity_median
    hypersensitivity_CV
    queen_allele_mean
    queen_allele_median
    queen_allele_CV
    nurse_allele_mean
    nurse_allele_median
    nurse_allele_CV
    forager_allele_mean
    forager_allele_median
    forager_allele_CV
    healthy_queens_number
    infected_queens_number
    healthy_nurses_number
    infected_nurses_number
    healthy_foragers_number
    infected_foragers_number
    infection_percentage
    total_infected
    total_bees
    died_by_aging
    died_by_non_aging
    died_by_pathogen
    died_by_hypersensitivity

    Records statistics for one simulation day and appends to the CSV file.
    Call this function only after run_id has been created.
    """
    try:    
        create_output_directory()
        var output_dir = get_output_directory()
        
        var os_module = Python.import_module("os")
        var builtins = Python.import_module("builtins")
        
        var csv_filename = simulation_data.run_id + "_dayStatistics.csv"
        var csv_path = os_module.path.join(output_dir, csv_filename)

        var day = simulation_data.days - 1
        if day < 0 or day >= len(simulation_data.daily_allele_data):
            print("Warning: No allele data for day", day)
            return
            
        if day >= len(simulation_data.daily_data):
            print("Warning: No population data for day", day)
            return

        var day_data = simulation_data.daily_allele_data[day].copy()

        # Get allele and fitness data from the day
        var fitness_values = List[Float64]()
        var queen_allele_values = List[Float64]()
        var nurse_allele_values = List[Float64]()
        var forager_allele_values = List[Float64]()
        var hypersensitivity_allele_values = List[Float64]()

        for swarm_data in day_data:
            fitness_values.append(Float64(swarm_data.fitness))
            queen_allele_values.append(Float64(swarm_data.queen_age_allele))
            nurse_allele_values.append(Float64(swarm_data.transition_age_allele))
            forager_allele_values.append(Float64(swarm_data.forager_age_allele))
            hypersensitivity_allele_values.append(Float64(swarm_data.hypersensitivity_allele))

        # Assign statistics storages
        var fitness_mean = 0.0
        var fitness_median: Float64
        var fitness_CV: Float64

        # Assign statistics storages
        var hypersensitivity_allele_mean = 0.0
        var hypersensitivity_allele_median: Float64
        var hypersensitivity_allele_CV: Float64

        var queen_allele_mean = 0.0
        var queen_allele_median: Float64
        var queen_allele_CV: Float64
        
        var nurse_allele_mean = 0.0
        var nurse_allele_median: Float64
        var nurse_allele_CV: Float64

        var forager_allele_mean = 0.0
        var forager_allele_median: Float64
        var forager_allele_CV: Float64

        var healthy_queens_number: Int
        var infected_queens_number: Int

        var healthy_nurses_number: Int
        var infected_nurses_number: Int

        var healthy_foragers_number: Int
        var infected_foragers_number: Int

        var infection_percentage: Float64
        var total_infected: Int
        var total_bees: Int

        var died_by_aging = 0
        var died_by_non_aging = 0
        var died_by_pathogen = 0
        var died_by_hypersensitivity = 0
        
        ### Calculate Metrics ###
        for fitness in fitness_values:
            fitness_mean += fitness
        fitness_mean /= Float64(len(fitness_values))
        fitness_median = calculate_median(fitness_values)

        for allele in queen_allele_values:
            queen_allele_mean += allele
        queen_allele_mean /= Float64(len(queen_allele_values))
        queen_allele_median = calculate_median(queen_allele_values)
        
        for allele in nurse_allele_values:
            nurse_allele_mean += allele
        nurse_allele_mean /= Float64(len(nurse_allele_values))
        nurse_allele_median = calculate_median(nurse_allele_values)
        
        for allele in forager_allele_values:
            forager_allele_mean += allele
        forager_allele_mean /= Float64(len(forager_allele_values))
        forager_allele_median = calculate_median(forager_allele_values)

        for allele in hypersensitivity_allele_values:
            hypersensitivity_allele_mean += allele
        hypersensitivity_allele_mean /= Float64(len(hypersensitivity_allele_values))
        hypersensitivity_allele_median = calculate_median(hypersensitivity_allele_values)
        
        var diff_sum = 0.0
        # Calculate coefficient variation CV
        for fitness in fitness_values:
            var diff = fitness - fitness_mean
            diff_sum += diff * diff
        var fitness_sd = sqrt(diff_sum / (Float64(len(fitness_values)) - 1) ) 
        fitness_CV = fitness_sd / fitness_mean

        for value in queen_allele_values:
            var diff = value - queen_allele_mean
            diff_sum += diff * diff
        var queen_allele_sd = sqrt(diff_sum / (Float64(len(queen_allele_values)) - 1) ) 
        queen_allele_CV = queen_allele_sd / queen_allele_mean

        for value in nurse_allele_values:
            var diff = value - nurse_allele_mean
            diff_sum += diff * diff
        var nurse_allele_sd = sqrt(diff_sum / (Float64(len(nurse_allele_values)) - 1) ) 
        nurse_allele_CV = nurse_allele_sd / nurse_allele_mean

        for value in forager_allele_values:
            var diff = value - forager_allele_mean
            diff_sum += diff * diff
        var forager_allele_sd = sqrt(diff_sum / (Float64(len(forager_allele_values)) - 1) ) 
        forager_allele_CV = forager_allele_sd / forager_allele_mean

        for value in hypersensitivity_allele_values:
            var diff = value - hypersensitivity_allele_mean
            diff_sum += diff * diff
        var hypersensitivity_allele_sd = sqrt(diff_sum / (Float64(len(hypersensitivity_allele_values)) - 1) ) 
        hypersensitivity_allele_CV = hypersensitivity_allele_sd / hypersensitivity_allele_mean

        # Calculate healthy castes
        healthy_queens_number, healthy_nurses_number, healthy_foragers_number = count_healthy_castes(swarms) 

        # Calculate infected castes
        infected_queens_number, infected_nurses_number, infected_foragers_number = count_infected_castes(swarms) 

        # Calculate infection_proportion
        infection_percentage, total_infected, total_bees = calculate_final_infection_proportion(simulation_data)

        # Calculate how many bees died in all 300 swarms this day by different reasons
        
        # by age
        for swarm_aging_death_at_day in simulation_data.daily_aging_dead[day]:
            died_by_aging += swarm_aging_death_at_day

        # by non aging random 
        for swarm_non_aging_death_at_day in simulation_data.daily_non_aging_dead[day]:
            died_by_non_aging += swarm_non_aging_death_at_day

        # by pathogen
        for swarm_pathogen_death_at_day in simulation_data.daily_pathogen_dead[day]:
            died_by_pathogen += swarm_pathogen_death_at_day

        # by hypersensitivity
        for swarm_hypersensitivity_death_at_day in simulation_data.daily_hypersensitivity_dead[day]:
            died_by_hypersensitivity += swarm_hypersensitivity_death_at_day

        # Prepare header string (base headers)
        var base_headers = [
            "day",
            "fitness_mean",
            "fitness_median",
            "fitness_CV",
            "hypersensitivity_allele_mean",
            "hypersensitivity_allele_median",
            "hypersensitivity_allele_CV",
            "queen_allele_mean",
            "queen_allele_median",
            "queen_allele_CV",
            "nurse_allele_mean",
            "nurse_allele_median",
            "nurse_allele_CV",
            "forager_allele_mean",
            "forager_allele_median",
            "forager_allele_CV",
            "healthy_queens_number",
            "infected_queens_number",
            "healthy_nurses_number",
            "infected_nurses_number",
            "healthy_foragers_number",
            "infected_foragers_number",
        ]

        # Build header string (include pathogen-specific columns)
        var header_string = ",".join(base_headers)
        for pathogen_idx in range(len(pathogens)):
            header_string += "," + "infected_by_pathogen" + String(pathogen_idx)
            header_string += "," + "immune_by_pathogen" + String(pathogen_idx)

        header_string += ",infection_percentage,total_infected,total_bees,died_by_aging,died_by_non_aging,died_by_pathogen,died_by_hypersensitivity\n"

        # Prepare row values in the same order as header
        var values_to_add = List[String]()
        values_to_add.append(String(day))
        values_to_add.append(String(fitness_mean))
        values_to_add.append(String(fitness_median))
        values_to_add.append(String(fitness_CV))
        values_to_add.append(String(hypersensitivity_allele_mean))
        values_to_add.append(String(hypersensitivity_allele_median))
        values_to_add.append(String(hypersensitivity_allele_CV))
        values_to_add.append(String(queen_allele_mean))
        values_to_add.append(String(queen_allele_median))
        values_to_add.append(String(queen_allele_CV))
        values_to_add.append(String(nurse_allele_mean))
        values_to_add.append(String(nurse_allele_median))
        values_to_add.append(String(nurse_allele_CV))
        values_to_add.append(String(forager_allele_mean))
        values_to_add.append(String(forager_allele_median))
        values_to_add.append(String(forager_allele_CV))
        values_to_add.append(String(healthy_queens_number))
        values_to_add.append(String(infected_queens_number))
        values_to_add.append(String(healthy_nurses_number))
        values_to_add.append(String(infected_nurses_number))
        values_to_add.append(String(healthy_foragers_number))
        values_to_add.append(String(infected_foragers_number))

        # Append per-pathogen values from simulation_data (safe access)
        var infected_counts = List[Int]()
        var immune_counts = List[Int]()
        if day < len(simulation_data.infected_by_pathogen):
            infected_counts = simulation_data.infected_by_pathogen[day].copy()
        if day < len(simulation_data.immune_to_pathogen):
            immune_counts = simulation_data.immune_to_pathogen[day].copy()

        for pathogen_idx in range(len(pathogens)):
            var inf_val = 0
            var imm_val = 0
            if pathogen_idx < len(infected_counts):
                inf_val = infected_counts[pathogen_idx]
            if pathogen_idx < len(immune_counts):
                imm_val = immune_counts[pathogen_idx]
            values_to_add.append(String(inf_val))
            values_to_add.append(String(imm_val))

        # Append final metrics
        values_to_add.append(String(infection_percentage))
        values_to_add.append(String(total_infected))
        values_to_add.append(String(total_bees))
        values_to_add.append(String(died_by_aging))
        values_to_add.append(String(died_by_non_aging))
        values_to_add.append(String(died_by_pathogen))
        values_to_add.append(String(died_by_hypersensitivity))

        # Create CSV line from values
        var csv_line = String("")
        for i in range(len(values_to_add)):
            csv_line += values_to_add[i]
            if i < len(values_to_add) - 1:
                csv_line += ","
        csv_line += "\n"
        
        # Check if file exists to write header
        var file_exists = os_module.path.exists(csv_path)
        
        var file_mode = "a"  # Append mode
        if not file_exists:
            file_mode = "w"  # Write mode (creates file)
        
        var file = builtins.open(csv_path, file_mode)
        
        # Write header if new file
        # Write header if new file
        if not file_exists:
            # header_string is built above to include per-pathogen headers
            _ = file.write(header_string)

        # Write data row
        _ = file.write(csv_line)
        _ = file.close()

    except e:
        pass


### For heatmap analysis ###
fn collect_beta_results(simulation_data: SimulationData, swarms: List[Swarm], pathogens: List[Pathogen], beta: Float64):
    """Collect results for beta analysis for one simulation and append to CSV file."""
    try:
        create_output_directory()
        var output_dir = get_output_directory()
        
        var os_module = Python.import_module("os")
        var builtins = Python.import_module("builtins")
        
        # Extract beta value from run_id or use provided beta
        var actual_beta = beta
        if actual_beta == 0.0 or len(pathogens) == 0:
            actual_beta = 0.0  # No pathogens case
        
        # Create filename based on beta value
        var beta_str = String(actual_beta)
        var csv_filename = "beta_" + beta_str + "_results.csv"
        var csv_path = os_module.path.join(output_dir, csv_filename)
        
        print("Collecting beta results for file:", csv_path)
        
        # Extract allele distribution from run_id - FIXED VERSION
        var distribution_info: String = ""
        var run_id_str = String(simulation_data.run_id)
        
        print("Run ID:", run_id_str)
        
        # Parse the distribution numbers from run_id
        # Format example: "pathogen_number_0_W_test_techSample1_dist1_1000-1000_32-32_10-100"
        if "_dist" in run_id_str:
            var parts = run_id_str.split("_")
            
            # Find the part that starts with numbers (queen_min-queen_max)
            for i in range(len(parts)):
                var part_str = String(parts[i])
                
                # Look for parts that contain "-" (like "1000-1000", "32-32", "10-100")
                if "-" in part_str:
                    # Check if this is likely a distribution part (contains only numbers and dash)
                    var is_distribution_part = True
                    for char in part_str.codepoints():
                        if String(char) != "-" and (String(char) < "0" or String(char) > "9"):
                            is_distribution_part = False
                            break
                    
                    if is_distribution_part:
                        # Found a distribution part, now collect all four parts
                        var distribution_parts = List[String]()
                        
                        # Get current and next two parts that match the pattern
                        if i < len(parts) and "-" in String(parts[i]):
                            distribution_parts.append(String(parts[i]))
                        if i + 1 < len(parts) and "-" in String(parts[i + 1]):
                            distribution_parts.append(String(parts[i + 1]))
                        if i + 2 < len(parts) and "-" in String(parts[i + 2]):
                            distribution_parts.append(String(parts[i + 2]))
                        if i + 3 < len(parts) and "-" in String(parts[i + 3]):
                            distribution_parts.append(String(parts[i + 3]))
                        
                        # If we found exactly 4 distribution parts, combine them
                        if len(distribution_parts) == 4:
                            distribution_info = distribution_parts[0] + "_" + distribution_parts[1] + "_" + distribution_parts[2] + "_" + distribution_parts[3]
                            print("Found distribution:", distribution_info)
                            break
        
        # Calculate statistics from the last day
        if simulation_data.days == 0:
            print("No simulation data available")
            return
        
        var last_day = simulation_data.days - 1
        
        # Get allele data from last day
        var fitness_values = List[Float64]()
        var nurse_allele_values = List[Float64]()
        var forager_allele_values = List[Float64]()
        var hypersensitivity_allele_values = List[Float64]()
        
        if last_day < len(simulation_data.daily_allele_data):
            var last_day_data = simulation_data.daily_allele_data[last_day].copy()
            
            for swarm_data in last_day_data:
                fitness_values.append(Float64(swarm_data.fitness))
                nurse_allele_values.append(Float64(swarm_data.transition_age_allele))
                forager_allele_values.append(Float64(swarm_data.forager_age_allele))
                hypersensitivity_allele_values.append(Float64(swarm_data.hypersensitivity_allele))
                
        # Calculate statistics
        var fitness_avg = 0.0
        var fitness_CV = 0.0
        var fitness_median = 0.0

        var nurse_allele_avg = 0.0
        var forager_allele_avg = 0.0
        var hypersensitivity_allele_avg = 0.0
        
        var infection_percentage = 0.0
        
        if len(fitness_values) > 0:
            # Calculate fitness mean, CV and median
            for fitness in fitness_values:
                fitness_avg += fitness
            fitness_avg /= Float64(len(fitness_values))
            fitness_median = calculate_median(fitness_values)
            
            for allele in nurse_allele_values:
                nurse_allele_avg += allele
            nurse_allele_avg /= Float64(len(nurse_allele_values))
            
            for allele in forager_allele_values:
                forager_allele_avg += allele
            forager_allele_avg /= Float64(len(forager_allele_values))

            for allele in hypersensitivity_allele_values:
                hypersensitivity_allele_avg += allele
            hypersensitivity_allele_avg /= Float64(len(hypersensitivity_allele_values))
            
            var diff_sum = 0.0
            # Calculate coefficient variation CV
            for fitness in fitness_values:
                var diff = fitness - fitness_avg
                diff_sum += diff * diff
            var fitness_sd = sqrt(diff_sum / (Float64(len(fitness_values)) - 1) ) 
            fitness_CV = fitness_sd / fitness_avg

            #count_infected_castes()
            # Calculate infection_proportion
            infection_percentage = calculate_final_infection_proportion(simulation_data)[0]

        # Prepare data row
        var row_data = List[String]()
        row_data.append(distribution_info)                    # Distribution identifier
        row_data.append(String(fitness_avg))                 # 1. Fitness average
        row_data.append(String(fitness_CV))            # 2. Fitness coefficient variation 
        row_data.append(String(fitness_median))       # Median
        row_data.append(String(nurse_allele_avg))            # 3. Nurse allele value
        row_data.append(String(forager_allele_avg))          # 4. Forager allele value
        row_data.append(String(hypersensitivity_allele_avg))
        row_data.append(String(infection_percentage))  # Infection percentage
        
        # Create CSV line
        var csv_line = String("")
        for i in range(len(row_data)):
            csv_line += row_data[i]
            if i < len(row_data) - 1:
                csv_line += ","
        csv_line += "\n"
        
        # Check if file exists to write header
        var file_exists = os_module.path.exists(csv_path)
        
        var file_mode = "a"  # Append mode
        if not file_exists:
            file_mode = "w"  # Write mode (creates file)
        
        var file = builtins.open(csv_path, file_mode)
        
        # Write header if new file
        if not file_exists:
            _ = file.write("distribution,fitness_avg,fitness_CV,fitness_median,nurse_allele_avg,forager_allele_avg,hypersensitivity_allele_avg,infection_percentage\n")
        
        # Write data row
        _ = file.write(csv_line)
        _ = file.close()
        
        print("✓ Beta results appended to:", csv_path)
        print("  Distribution:", distribution_info)
        print("  Fitness avg:", fitness_avg, "coefficient variation:", fitness_CV)
        print("  Fitness median:", fitness_median, "infection percentage:", infection_percentage)
        print("  Nurse allele:", nurse_allele_avg, "Forager allele:", forager_allele_avg)
        print("  Hypersensitivity allele:", hypersensitivity_allele_avg)
        
    except e:
        print("Failed to collect beta results. Error:", e)








