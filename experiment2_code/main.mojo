from types import *
from simulation_core import *
from pathogen_system import *
from analysis import *
from randomizer import Xoshiro256PlusPlus
from time import perf_counter_ns
from random import seed
from math import ceil
from algorithm.functional import parallelize

# Set parameters before running
comptime MUTATION_RATE: Float64 = 0.1 # influences mutation chance in function inside simulation_core.mojo
comptime RANDOM_SWARM_DEATH: Float64 = 0.0005
comptime MAX_SIMULATION_DAYS: Int = 356000
comptime FORAGER_DEATH_RATE: Float64 = 0.134

comptime NUM_REPLICATES: Int = 10 # Number of times to repeat each experiment, minimum 1


fn main():
    """Main simulation function."""
    seed(42) # set seed to update randomizator

    # Define different allele distributions to test - NEW FIXED VALUES
    var allele_distributions = List[Tuple[Int, Int, Int, Int, Int, Int, Int, Int]]()
    
    allele_distributions.append((1, 1000, 1, 60, 1, 1000, 0, 50))

    var num_of_pathogens = 1

    # Define beta values conditionally
    var beta_values = List[Float64]()
    if num_of_pathogens == 0:
        beta_values.append(0.0)     # dummy, ensures loop executes once
    else:
        # beta_values.append(0.001)
        # beta_values.append(0.01)
        # beta_values.append(0.05)
        # beta_values.append(0.08)
        # beta_values.append(0.2)
        # beta_values.append(0.19)
        # beta_values.append(0.13)
        # beta_values.append(0.0)
        # beta_values.append(0.27)
        # beta_values.append(0.3)
        # beta_values.append(0.4)
        # beta_values.append(0.5)
        # beta_values.append(0.9)
        # beta_values.append(0.17)

        # beta_values.append(0.0)

        for b in range(0.0, 0.151, 0.01):
            beta_values.append(round(b, 3))

        # beta_values.append(0.1)


        # for b in range(0.0, 4.51, 0.225):
        #     beta_values.append(round(b, 3))

        # beta_values.append(4.5)


    # START FROM BEGINNING - NEW EXPERIMENT SET
    var resume_dist_idx = 0      # Start from dist1
    var resume_beta_idx = 0      # Start from beta 0.1
    var resume_replicate = 1  # Start from techSample1

    var start = perf_counter_ns()
    
    # Loop over each combination - STARTING FROM BEGINNING
    for dist_idx in range(resume_dist_idx, len(allele_distributions)):
        var distribution = allele_distributions[dist_idx]
        var queen_min = distribution[0]
        var queen_max = distribution[1] 
        var transition_min = distribution[2]
        var transition_max = distribution[3]
        var forager_min = distribution[4]
        var forager_max = distribution[5]
        var hypersensitivity_min = distribution[6]
        var hypersensitivity_max = distribution[7]
        
        
        # Outer loop over tech samples (replicates), inner loop over beta
        for replicate in range(resume_replicate, NUM_REPLICATES + 1):

            # Determine starting beta index for this replicate
            var start_beta_idx = 0
            if dist_idx == resume_dist_idx and replicate == resume_replicate:
                start_beta_idx = resume_beta_idx

            for beta_idx in range(start_beta_idx, len(beta_values)):
                var beta = beta_values[beta_idx]

                # Per-sample base seed incorporated into per-swarm RNGs below

                # Create run identifier with techSample prefix
                # Build run identifier conditionally on pathogen count
                var run_id = "pathogen_number_" + String(num_of_pathogens) +
                            "_W_test_techSample" + String(replicate) +
                            "_dist" + String(dist_idx + 1) + "_" +
                            String(queen_min) + "-" + String(queen_max) + "_" +
                            String(transition_min) + "-" + String(transition_max) + "_" +
                            String(forager_min) + "-" + String(forager_max) + "_" +
                            String(hypersensitivity_min) + "-" + String(hypersensitivity_max)

                # Append beta only when num_of_pathogens > 0
                if num_of_pathogens > 0:
                    run_id += "_beta" + String(beta)

                print("=" * 80)
                print("STARTING RUN: " + run_id)
                print("Distribution", dist_idx + 1, "- Queen:", queen_min, "-", queen_max, 
                      ", Transition:", transition_min, "-", transition_max,
                      ", Forager:", forager_min, "-", forager_max,
                      ", Hypersensitivity:", hypersensitivity_min, "-", hypersensitivity_max)
                print("Beta:", beta, "| Replicate:", replicate, "/", NUM_REPLICATES)
                
                # Create pathogens for this run
                var pathogens = List[Pathogen]()
                
                pathogens.append(Pathogen(0, beta, True, True, False, 0.0, 0.0))
                # pathogens.append(Pathogen(1, 4.5, False, False, False, 0.0, 0.0))


                # Create swarms with specific chromosome ranges
                var swarms = create_swarms_with_ranges(
                queen_min, queen_max,
                transition_min, transition_max,
                forager_min, forager_max,
                0, 0)

                # #set  which allele will be binary distributed (so only two values would be in the whole population)
                for i in range(NUM_SWARMS):
                    if i % 2 == 0:
                        swarms[i].chromosome.allele1[3] = hypersensitivity_min
                        swarms[i].chromosome.allele2[3] = hypersensitivity_min
                    else:
                        swarms[i].chromosome.allele1[3] = hypersensitivity_max
                        swarms[i].chromosome.allele2[3] = hypersensitivity_max
                    
                # Initialize data collection with run ID
                var simulation_data = SimulationData(run_id)
                
                # Record initial data
                simulation_data.add_day_data(swarms, pathogens)
                
                # Simulate aging until steady state is reached
                var day = 1
                
                # Initialize per-swarm RNGs (thread-local)
                var rngs = List[Xoshiro256PlusPlus]()
                for i in range(len(swarms)):
                    rngs.append(Xoshiro256PlusPlus(UInt64(42 + replicate * 9973 + i * 37)))
                
                print("Starting simulation until steady state is reached...")
                print("Maximum simulation days:", MAX_SIMULATION_DAYS)
                print()

                # store times of function execution
                var time_check_swarm_death = 0.0
                var time_process_eggs_and_hatching = 0.0
                var time_infect_foragers = 0.0
                var time_spread_infection = 0.0
                var time_age_and_kill_bees = 0.0

                var time_collect_day_statistics = 0.0
                var time_add_day_data = 0.0
                var time_check_steady_state = 0.0

                while day <= MAX_SIMULATION_DAYS:

                    # Prepare per-day death stats buffers and ensure day lists exist
                    var current_day = simulation_data.days - 1
                    if len(simulation_data.daily_aging_dead) <= current_day:
                        simulation_data.daily_aging_dead.append([])
                        simulation_data.daily_non_aging_dead.append([])
                        simulation_data.daily_pathogen_dead.append([])
                        simulation_data.daily_hypersensitivity_dead.append([])

                    var num_swarms_today = len(swarms)
                    var aging_counts = List[Int]()
                    var non_aging_counts = List[Int]()
                    var pathogen_counts = List[Int]()
                    var hypersensitivity_counts = List[Int]()
                    var t_check_swarm_death = List[Float64]()
                    var t_process_eggs = List[Float64]()
                    var t_infect_foragers = List[Float64]()
                    var t_spread_infection = List[Float64]()
                    var t_age_and_kill = List[Float64]()

                    for _ in range(num_swarms_today):
                        aging_counts.append(0)
                        non_aging_counts.append(0)
                        pathogen_counts.append(0)
                        hypersensitivity_counts.append(0)
                        t_check_swarm_death.append(0.0)
                        t_process_eggs.append(0.0)
                        t_infect_foragers.append(0.0)
                        t_spread_infection.append(0.0)
                        t_age_and_kill.append(0.0)

                    @parameter
                    fn process_swarm(i: Int):
                        ref swarm = swarms[i]

                        var start_csd = perf_counter_ns()
                        check_swarm_death(swarm, swarm_death_rate=RANDOM_SWARM_DEATH)
                        var end_csd = perf_counter_ns()
                        t_check_swarm_death[i] = count_time_sec(start_csd, end_csd)

                        if len(swarm.bees) > 0:
                            var start_eggs = perf_counter_ns()
                            process_eggs_and_hatching(swarm, pathogens)
                            var end_eggs = perf_counter_ns()
                            t_process_eggs[i] = count_time_sec(start_eggs, end_eggs)

                            if num_of_pathogens > 0:
                                var start_inf = perf_counter_ns()
                                infect_foragers(swarm, pathogens, rngs[i])
                                var end_inf = perf_counter_ns()
                                t_infect_foragers[i] = count_time_sec(start_inf, end_inf)

                                var start_spread = perf_counter_ns()
                                spread_infection(swarm, pathogens, rngs[i])
                                var end_spread = perf_counter_ns()
                                t_spread_infection[i] = count_time_sec(start_spread, end_spread)

                            var start_age = perf_counter_ns()
                            var stats = age_and_kill_bees(swarm, pathogens, simulation_data, rngs[i])
                            var end_age = perf_counter_ns()
                            t_age_and_kill[i] = count_time_sec(start_age, end_age)

                            aging_counts[i] = stats[0]
                            non_aging_counts[i] = stats[1]
                            pathogen_counts[i] = stats[2]
                            hypersensitivity_counts[i] = stats[3]

                    parallelize[process_swarm](num_swarms_today)

                    # Aggregate timings using MAX (wall time, not CPU time)
                    # For parallel sections, use the longest swarm's execution time
                    time_check_swarm_death += max_f64_list(t_check_swarm_death)
                    time_process_eggs_and_hatching += max_f64_list(t_process_eggs)
                    time_infect_foragers += max_f64_list(t_infect_foragers)
                    time_spread_infection += max_f64_list(t_spread_infection)
                    time_age_and_kill_bees += max_f64_list(t_age_and_kill)

                    # Append per-swarm death counts to simulation_data for current day
                    for i in range(num_swarms_today):
                        simulation_data.daily_aging_dead[current_day].append(aging_counts[i])
                        simulation_data.daily_non_aging_dead[current_day].append(non_aging_counts[i])
                        simulation_data.daily_pathogen_dead[current_day].append(pathogen_counts[i])
                        simulation_data.daily_hypersensitivity_dead[current_day].append(hypersensitivity_counts[i])

                    print('Day: ' + String(day))

                    # Record daily data (include pathogens so per-pathogen columns can be written)
                    var start_collect_day_statistics = perf_counter_ns()
                    if day % 10 == 0:
                        collect_day_statistics(simulation_data, swarms, pathogens)
                    var end_collect_day_statistics = perf_counter_ns()
                    time_collect_day_statistics += count_time_sec(start_collect_day_statistics, end_collect_day_statistics)

                    var start_add_day_data = perf_counter_ns()
                    simulation_data.add_day_data(swarms, pathogens)
                    var end_add_day_data = perf_counter_ns()
                    time_add_day_data += count_time_sec(start_add_day_data, end_add_day_data)
                    
                    # Check steady state
                    var start_time_check_steady_state = perf_counter_ns()
                    if check_steady_state(simulation_data,
                    coefficient_variation_threshold=0.1,
                    mean_difference_threshold=0.1,
                    lookback_days=100, min_stable_checks=5, day_check_window=100):
                        print()
                        print("🎯 SIMULATION COMPLETED: Steady state reached after", day, "days!")
                        break

                    var end_time_check_steady_state = perf_counter_ns()
                    time_check_steady_state  += count_time_sec(start_time_check_steady_state , end_time_check_steady_state)

                    # Handle reproduction for dead swarms
                    var reproductions_today = 0
                    for i in range(len(swarms)):
                        if len(swarms[i].bees) == 0:  # Dead swarm found
                            record_swarm_death(simulation_data)
                            var parent_indices = weighted_select_parents(swarms, pathogens)
                            
                            if parent_indices[0] != -1 and parent_indices[1] != -1:
                                var parent1_chromosome = swarms[parent_indices[0]].chromosome.copy()
                                var parent2_chromosome = swarms[parent_indices[1]].chromosome.copy()
                                reproduce_dead_swarm(swarms[i], parent1_chromosome, parent2_chromosome, rng=rngs[i])
                                reproductions_today += 1

                    day += 1

                if day > MAX_SIMULATION_DAYS:
                    print()
                    print("⚠️  SIMULATION STOPPED: Maximum simulation time reached")

                # Print final results for this run
                print()
                print("FINAL RESULTS FOR RUN:", run_id)
                print("Days simulated:", day - 1)
                
                # Create visualization for this run
                print()
                print("Creating visualizations for run:", run_id)
                # create_animated_visualization(simulation_data)
                print("✓ Visualizations completed for run:", run_id)

                # Collect beta results for analysis
                # collect_beta_results(simulation_data, swarms, pathogens, beta)
                
                print()
                print("RUN", run_id, "COMPLETED!")
                print()

                print("=" * 80)
                print("PERFORMANCE TIMING FOR THIS RUN SUMMARY")
                print("=" * 80)

                var total_function_time = (
                    time_check_swarm_death + 
                    time_process_eggs_and_hatching + 
                    time_infect_foragers + 
                    time_spread_infection + 
                    time_age_and_kill_bees + 
                    time_collect_day_statistics + 
                    time_add_day_data + 
                    time_check_steady_state
                )
                print('General function usage time: ' + String(total_function_time))
                print()
                print("check_swarm_death, seconds: " + String(time_check_swarm_death) + ", percent of general time: " + String(time_check_swarm_death / (total_function_time / 100)))
                print("process_eggs_and_hatching, seconds: " + String(time_process_eggs_and_hatching) + ", percent of general time: " + String(time_process_eggs_and_hatching / (total_function_time / 100)))
                print("infect_foragers, seconds: " + String(time_infect_foragers) + ", percent of general time: " + String(time_infect_foragers / (total_function_time / 100)))
                print("spread_infection, seconds: " + String(time_spread_infection) + ", percent of general time: " + String(time_spread_infection / (total_function_time / 100)))
                print("age_and_kill_bees, seconds: " + String(time_age_and_kill_bees) + ", percent of general time: " + String(time_age_and_kill_bees / (total_function_time / 100)))
                print("collect_day_statistics, seconds: " + String(time_collect_day_statistics) + ", percent of general time: " + String(time_collect_day_statistics / (total_function_time / 100)))
                print("add_day_data, seconds: " + String(time_add_day_data) + ", percent of general time: " + String(time_add_day_data / (total_function_time / 100)))
                print("check_steady_state, seconds: " + String(time_check_steady_state) + ", percent of general time: " + String(time_check_steady_state / (total_function_time / 100)))

                print("=" * 80)

                print('Number of days / iterations: ' + String(day))

    var end = perf_counter_ns()
    var total_time = end - start            
    var minutes_total_time = total_time / 1e+9 / 60
    var seconds_total_time = total_time / 1e+9
    print("=" * 80)
    print("ALL PARAMETER SWEEP RUNS COMPLETED!")
    print("Execution time, min: " + String(minutes_total_time))
    print("Execution time, sec: " + String(seconds_total_time))
    print("=" * 80)



fn max_f64_list(xs: List[Float64]) -> Float64:
    """Return maximum value from list of Float64."""
    var m = 0.0
    for x in xs:
        if x > m:
            m = x
    return m
