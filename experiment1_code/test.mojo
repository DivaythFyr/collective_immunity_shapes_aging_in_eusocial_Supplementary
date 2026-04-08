from random import random_float64, random_si64
from time import perf_counter_ns
from randomizer import Xoshiro256PlusPlus
from math import *
from pathogen_system import calculate_infection_shots
from simulation_core import create_swarms_with_ranges, generate_equally_distributed_values, reproduce_dead_swarm
from analysis import calculate_fitness
from analysis import calculate_median
from analysis import bubble_sort

from types import *
from simulation_core import attempt_cure
from pathogen_system import spread_infection

from simulation_core import create_swarms_with_ranges


fn main():
    # test_XOR()

    # divide_()

    # test_initial_distribution()

    # test_reproduce_mutations()


    test_mutation_distribution()


fn test_XOR():
    var rng = Xoshiro256PlusPlus(1234500)
    var StartTime = perf_counter_ns()

    var number: Float64
    for _ in range(1000000):
        number = rng.next_float()
        # number = rng.next_range(0, 1000000)
        if number < 0.0001:
            print(number)
        # print(number)
    var EndTime = perf_counter_ns()
    print("Time taken:", (EndTime - StartTime) / 1_000_000, "seconds")

fn divide_():
    var a = UInt32.MAX 
    var b = UInt32.MAX - 100000 
    # var c = Float64(a) / Float64(b)
    # var c = Float64((UInt128(a) * 1000000) // UInt128(b))# / 1000000
    var c = Float64((Float64(a) * 1000000) // Float64(b))# / 1000000
    print(c, 'a=' + String(a) + ' b=' + String(b))


    a = UInt32.MAX - 100000
    b = UInt32.MAX - 200000
    # c = Float64((UInt128(a) * 1000000) // UInt128(b))# / 1000000
    c = Float64((Float64(a) * 1000000) // Float64(b))# / 1000000
    # c = Float64(a) / Float64(b)
    
    print(c, 'a=' + String(a) + ' b=' + String(b))

    var f1 = 1.0000000001
    var f2 = 1.0000000002
    var f3 = f1 + f2

    print(f3, 'f1=' + String(f1) + ' f2=' + String(f2))


    print(UInt32.MAX)
    print(Float64(UInt32.MAX))


    f1 = 0.000001
    f2 = 0.000001
    f3 = f1 * f2
    print(f1, ' * ',  f2,  ' = ', f3)


    f1 = 1.0 
    f2 = 4294967295.0
    f3 = f1 * f2
    print(f1, ' * ',  f2,  ' = ', f3)


    f1 = 1.0 
    f2 = 4294967295.0
    f3 = f1 / f2
    print(f1, ' / ',  f2,  ' = ', f3)

    f1 = 4294967295.0 
    f2 = 4294967294.0
    f3 = f1 / f2
    print(f1, ' / ',  f2,  ' = ', f3)


    f1 = 4294967294.0 
    f2 = 4294967293.0
    f3 = f1 / f2
    print(f1, ' / ',  f2,  ' = ', f3)


fn test_initial_distribution():
    # Distribution: 1-1000_1-60_1-1000_0-50
    # Gene 1 (Eusociality allele):      min=1,   max=1000
    # Gene 2 (Lifespan allele):         min=1,   max=60
    # Gene 3 (Immunity allele):         min=1,   max=1000
    # Gene 4 (Hypersensitivity allele): min=0,   max=50

    var queen_allele_min: Int = 1
    var queen_allele_max: Int = 1000
    var nurse_allele_min: Int = 1
    var nurse_allele_max: Int = 60
    var forager_allele_min: Int = 1
    var forager_allele_max: Int = 1000
    var hypersensitivity_min: Int = 0
    var hypersensitivity_max: Int = 50

    print("Initializing swarms with distribution: 1-1000_1-60_1-1000_0-50")
    print("  Queen allele range:      " + String(queen_allele_min) + " - " + String(queen_allele_max))
    print("  Nurse allele range:         " + String(nurse_allele_min) + " - " + String(nurse_allele_max))
    print("  Forager allele range:         " + String(forager_allele_min) + " - " + String(forager_allele_max))
    print("  Hypersensitivity allele range: " + String(hypersensitivity_min) + " - " + String(hypersensitivity_max))

    var pathogens = List[Pathogen]()
                
    pathogens.append(Pathogen(0, 0.1, True, True, False, 0.0, 0.0))
    pathogens.append(Pathogen(1, 4.5, False, False, False, 0.0, 0.0))

    var swarms = create_swarms_with_ranges(
        queen_allele_min, queen_allele_max,
        nurse_allele_min, nurse_allele_max,
        forager_allele_min, forager_allele_max,
        hypersensitivity_min, hypersensitivity_max,
    )

    print("\nNumber of swarms created: " + String(len(swarms)))

    # Print allele distribution for each swarm
    print("\nAllele distribution across all swarms:")
    print("Swarm | Queen | Nurse | Forager | Hypersensitivity")
    print("------+-------------+----------+----------+-----------------")
    for i in range(len(swarms)):
        var eusociality_val = swarms[i].chromosome.get_trait_value(0)
        var lifespan_val = swarms[i].chromosome.get_trait_value(1)
        var immunity_val = swarms[i].chromosome.get_trait_value(2)
        var hypersensitivity_val = swarms[i].chromosome.get_trait_value(3)
        print(
            "  " + String(i) 
            + "   |     " + String(eusociality_val) 
            + "     |    " + String(lifespan_val) 
            + "    |    " + String(immunity_val) 
            + "    |       " + String(hypersensitivity_val)
        )


fn test_reproduce_mutations():
    """Test how often mutations occur during reproduce_dead_swarm 
    and what value changes they cause."""
    
    print("=" * 60)
    print("TEST: Mutation frequency and value changes in reproduce_dead_swarm")
    print("=" * 60)
    
    var num_trials = 1000000
    var mutation_count = 0
    
    # Track mutation details: gene index, old value, new value, delta
    var mutated_genes = List[Int]()
    var mutation_deltas = List[Int]()
    var mutation_old_values = List[Int]()
    var mutation_new_values = List[Int]()
    
    var rng = Xoshiro256PlusPlus(42)
    
    # Create two parent chromosomes with known values
    var parent1_chromosome = Chromosome()
    parent1_chromosome.allele1.append(500)   # Queen age gene
    parent1_chromosome.allele2.append(500)
    parent1_chromosome.allele1.append(30)    # Transition age gene
    parent1_chromosome.allele2.append(30)
    parent1_chromosome.allele1.append(500)   # Forager age gene
    parent1_chromosome.allele2.append(500)
    parent1_chromosome.allele1.append(25)    # Hypersensitivity gene
    parent1_chromosome.allele2.append(25)
    
    var parent2_chromosome = Chromosome()
    parent2_chromosome.allele1.append(600)   # Queen age gene
    parent2_chromosome.allele2.append(600)
    parent2_chromosome.allele1.append(40)    # Transition age gene
    parent2_chromosome.allele2.append(40)
    parent2_chromosome.allele1.append(600)   # Forager age gene
    parent2_chromosome.allele2.append(600)
    parent2_chromosome.allele1.append(30)    # Hypersensitivity gene
    parent2_chromosome.allele2.append(30)
    
    print("\nParent 1 chromosome: Queen=500, Transition=30, Forager=500, Hyper=25")
    print("Parent 2 chromosome: Queen=600, Transition=40, Forager=600, Hyper=30")
    print("\nRunning", num_trials, "reproduction trials...")
    print("Expected mutation rate (MUTATION_RATE):", MUTATION_RATE)
    print("")
    
    for trial in range(num_trials):
        # Create a dead swarm (empty bees list)
        var dead_swarm = Swarm()
        dead_swarm.chromosome = parent1_chromosome.copy()
        
        # Record offspring chromosome BEFORE mutation by doing crossover manually
        # We can't intercept, so instead compare offspring to expected parent values
        reproduce_dead_swarm(dead_swarm, parent1_chromosome, parent2_chromosome, rng)
        
        # Check each gene: if the value is NOT one of the parent alleles, mutation occurred
        for gene_idx in range(NUM_GENES):
            var offspring_allele1 = dead_swarm.chromosome.allele1[gene_idx]
            var offspring_allele2 = dead_swarm.chromosome.allele2[gene_idx]
            
            var p1_a1 = parent1_chromosome.allele1[gene_idx]
            var p1_a2 = parent1_chromosome.allele2[gene_idx]
            var p2_a1 = parent2_chromosome.allele1[gene_idx]
            var p2_a2 = parent2_chromosome.allele2[gene_idx]
            
            # Check allele1 - should come from parent1
            var possible_values_allele1 = List[Int]()
            possible_values_allele1.append(p1_a1)
            possible_values_allele1.append(p1_a2)
            
            var allele1_is_mutated = True
            for v in range(len(possible_values_allele1)):
                if offspring_allele1 == possible_values_allele1[v]:
                    allele1_is_mutated = False
                    break
            
            if allele1_is_mutated:
                mutation_count += 1
                mutated_genes.append(gene_idx)
                # Find closest parent value as "old value"
                var old_val = p1_a1
                var delta = offspring_allele1 - old_val
                mutation_old_values.append(old_val)
                mutation_new_values.append(offspring_allele1)
                mutation_deltas.append(delta)
            
            # Check allele2 - should come from parent2
            var possible_values_allele2 = List[Int]()
            possible_values_allele2.append(p2_a1)
            possible_values_allele2.append(p2_a2)
            
            var allele2_is_mutated = True
            for v in range(len(possible_values_allele2)):
                if offspring_allele2 == possible_values_allele2[v]:
                    allele2_is_mutated = False
                    break
            
            if allele2_is_mutated:
                mutation_count += 1
                mutated_genes.append(gene_idx)
                var old_val2 = p2_a1
                var delta2 = offspring_allele2 - old_val2
                mutation_old_values.append(old_val2)
                mutation_new_values.append(offspring_allele2)
                mutation_deltas.append(delta2)
    
    # Print results
    print("=" * 60)
    print("RESULTS")
    print("=" * 60)
    print("Total trials:          ", num_trials)
    print("Total mutations found: ", mutation_count)
    print("Mutation frequency:    ", Float64(mutation_count) / Float64(num_trials), "per reproduction")
    print("")
    
    # Count mutations per gene
    var gene_names = List[String]()
    gene_names.append("Queen Age")
    gene_names.append("Transition Age")
    gene_names.append("Forager Age")
    gene_names.append("Hypersensitivity")
    
    print("Mutations per gene:")
    for g in range(NUM_GENES):
        var count = 0
        for m in range(len(mutated_genes)):
            if mutated_genes[m] == g:
                count += 1
        print("  " + gene_names[g] + ": " + String(count) + " mutations")
    print("")
    
    # Print mutation delta statistics
    if len(mutation_deltas) > 0:
        var min_delta = mutation_deltas[0]
        var max_delta = mutation_deltas[0]
        var sum_delta: Float64 = 0.0
        var sum_abs_delta: Float64 = 0.0
        var positive_count = 0
        var negative_count = 0
        
        for i in range(len(mutation_deltas)):
            var d = mutation_deltas[i]
            sum_delta += Float64(d)
            if d < 0:
                sum_abs_delta += Float64(-d)
                negative_count += 1
            else:
                sum_abs_delta += Float64(d)
                positive_count += 1
            if d < min_delta:
                min_delta = d
            if d > max_delta:
                max_delta = d
        
        var avg_delta = sum_delta / Float64(len(mutation_deltas))
        var avg_abs_delta = sum_abs_delta / Float64(len(mutation_deltas))
        
        print("Mutation value change statistics:")
        print("  Min delta:              ", min_delta)
        print("  Max delta:              ", max_delta)
        print("  Average delta:          ", avg_delta)
        print("  Average absolute delta: ", avg_abs_delta)
        print("  Positive mutations:     ", positive_count)
        print("  Negative mutations:     ", negative_count)
        print("")
        
        # Print first 20 individual mutations for inspection
        var show_count = min(20, len(mutation_deltas))
        print("First", show_count, "mutations:")
        print("  Gene           | Old Value | New Value | Delta")
        print("  ---------------+-----------+-----------+------")
        for i in range(show_count):
            print(
                "  " + gene_names[mutated_genes[i]]
                + " | " + String(mutation_old_values[i])
                + " | " + String(mutation_new_values[i])
                + " | " + String(mutation_deltas[i])
            )
    else:
        print("No mutations detected in", num_trials, "trials.")
    
    print("\n" + "=" * 60)
    print("NOTE: mutate_chromosome excludes hypersensitivity gene (gene 3)")
    print("      Mutation range: ±1 to ±1000 per allele")
    print("=" * 60)


fn test_mutation_distribution():
    """Generate and print 100 mutation values using log10 distribution, then verify the distribution statistically."""
    var rng = Xoshiro256PlusPlus(42)  # Use a fixed seed for reproducibility
    var min_mutation_amount: Int = 1
    var max_mutation_amount: Int = 1000
    
    print("=" * 60)
    print("TEST: Log10 Mutation Distribution (100 samples)")
    print("=" * 60)
    print("Min mutation amount:", min_mutation_amount)
    print("Max mutation amount:", max_mutation_amount)
    print("Distribution: Smaller values are more likely (log10 bias)")
    print("")
    
    var log_min = log10(Float64(min_mutation_amount))
    var log_max = log10(Float64(max_mutation_amount))
    
    var log_values = List[Float64]()  # Collect log-transformed values for stats
    
    for i in range(100):
        var u = rng.next_float64()  # Uniform [0, 1)
        var log_val = log_min + (log_max - log_min) * u
        var raw_mutation_value = pow(10.0, log_val)
        var mutation_value = Int(round(raw_mutation_value))
        
        # Clamp to bounds
        mutation_value = max(min_mutation_amount, min(max_mutation_amount, mutation_value))
        
        # Collect log10 of the mutation value for statistical check
        log_values.append(log10(Float64(mutation_value)))
        
        print("Mutation", i + 1, ":", mutation_value)
    
    print("\n" + "=" * 60)
    print("Visual inspection: Expect more values near 1-10, fewer near 100-1000.")
    print("=" * 60)
    
    # Statistical check: Verify log-transformed values are uniform
    var expected_mean = (log_min + log_max) / 2.0
    var expected_variance = (log_max - log_min) ** 2 / 12.0
    
    var sum_log = 0.0
    var sum_sq_log = 0.0
    for lv in log_values:
        sum_log += lv
        sum_sq_log += lv * lv
    
    var n = Float64(len(log_values))
    var mean_log = sum_log / n
    var variance_log = (sum_sq_log / n) - (mean_log * mean_log)
    
    print("\nDistribution Verification (Log Scale):")
    print("Expected mean (uniform in log space):", expected_mean)
    print("Actual mean (log of mutations):      ", mean_log)
    print("Expected variance (uniform in log space):", expected_variance)
    print("Actual variance (log of mutations):      ", variance_log)
    print("If mean and variance are close to expected, the distribution is log-uniform.")
    print("=" * 60)

    












    



    









    

    
    