from types import *
from pathogen_system import is_queen_alive_and_fertile
from random import random_si64, shuffle
from main import *
from algorithm.functional import parallelize
from randomizer import Xoshiro256PlusPlus
from math import log10 

fn generate_base_longevity(caste_type: Int, chromosome: Chromosome) -> Int:
    """Generate base longevity based on caste type."""
    if caste_type == QUEEN:
        # Queens: use the chromosome-determined longevity
        result = chromosome.get_trait_value(QUEEN_AGE_GENE)
    elif caste_type == NURSE:
        # Nurses: transition age
        result = chromosome.get_trait_value(TRANSITION_AGE_GENE)
    else: # FORAGER
        # Foragers: use chromosome-determined age
        result = chromosome.get_trait_value(FORAGER_AGE_GENE)
    return result

fn generate_swarm_with_chromosome(chromosome: Chromosome) -> Swarm:
    """Generate a swarm with a specific chromosome.
     All bees will receive base logevity based on the chromosome.
     Individual initial age is set randomly from min to max."""
    var swarm = Swarm()
    swarm.chromosome = chromosome.copy()
    
    # Add queen(s)
    for _ in range(INITIAL_QUEENS):
        swarm.bees.append(Bee(caste_type = QUEEN,
         base_longevity = generate_base_longevity(QUEEN, swarm.chromosome),
         age=Int(random_si64(0, generate_base_longevity(QUEEN, swarm.chromosome)))
         ))
    # Add nurses
    for _ in range(INITIAL_NURSES):
        swarm.bees.append(Bee(caste_type = NURSE,
         base_longevity = generate_base_longevity(NURSE, swarm.chromosome),
         age=Int(random_si64(0, generate_base_longevity(NURSE, swarm.chromosome)))
         ))
    # Add foragers
    for _ in range(INITIAL_FORAGERS):
        swarm.bees.append(Bee(caste_type = FORAGER,
         base_longevity = generate_base_longevity(FORAGER, swarm.chromosome),
         age=Int(random_si64(0, generate_base_longevity(FORAGER, swarm.chromosome)))
         ))
    
    return swarm.copy()

fn create_swarms_with_ranges_random_distribution(queen_min: Int, queen_max: Int, transition_min: Int, transition_max: Int, forager_min: Int, forager_max: Int, hypersensitivity_min: Int, hypersensitivity_max: Int) -> List[Swarm]:
    """Create swarms with chromosomes, the allele values are given randomly from min to max."""
    var swarms = List[Swarm]()
    for _ in range(NUM_SWARMS):
        var chromosome = Chromosome()
        chromosome.generate_alleles(queen_min, queen_max, transition_min, transition_max, forager_min, forager_max, hypersensitivity_min, hypersensitivity_max)
        swarms.append(generate_swarm_with_chromosome(chromosome))
    return swarms.copy()    

fn generate_equally_distributed_values(min_val: Int, max_val: Int, num_samples: Int) -> List[Int]:
    """Generate num_samples equally distributed integer values between min_val and max_val. Used in swarm creation when total possible allele combinations is higher than number of swarms."""
    var values = List[Int]()
    
    if min_val == max_val:
        # Fixed value case - return the same value num_samples times
        for _ in range(num_samples):
            values.append(min_val)
    else:
        # Calculate step size and generate equally spaced values
        var range_size = max_val - min_val
        var step = Float64(range_size) / Float64(num_samples - 1)
        
        for i in range(num_samples):
            var float_val = Float64(min_val) + Float64(i) * step
            var int_val = Int(round(float_val))
            # Ensure we stay within bounds
            int_val = max(min_val, min(max_val, int_val))
            values.append(int_val)
    
    return values.copy()

fn create_swarms_with_ranges(queen_min: Int, queen_max: Int,
 transition_min: Int, transition_max: Int, forager_min: Int, forager_max: Int, hypersensitivity_min: Int, hypersensitivity_max: Int) -> List[Swarm]:
    """Create swarms with chromosomes ensuring equal distribution of allele values - NEW SAMPLING VERSION."""
    var swarms = List[Swarm]()
    
    # Calculate ranges and number of unique values for each allele
    var queen_range = queen_max - queen_min
    var transition_range = transition_max - transition_min  
    var forager_range = forager_max - forager_min
    var hypersensitivity_range = hypersensitivity_max - hypersensitivity_min
    
    var queen_unique_values = max(1, queen_range + 1)
    var transition_unique_values = max(1, transition_range + 1)
    var forager_unique_values = max(1, forager_range + 1)
    var hypersensitivity_unique_values = max(1, hypersensitivity_range + 1)
    
    print("Creating swarms with distribution:")
    print("  Queen alleles:", queen_unique_values, "values from", queen_min, "to", queen_max)
    print("  Nurse alleles:", transition_unique_values, "values from", transition_min, "to", transition_max)
    print("  Forager alleles:", forager_unique_values, "values from", forager_min, "to", forager_max)
    print("  Hypersensitivity alleles:", hypersensitivity_unique_values, "values from", hypersensitivity_min, "to", hypersensitivity_max)
    
    # Calculate how to distribute swarms
    var total_unique_combinations = queen_unique_values * transition_unique_values * forager_unique_values * hypersensitivity_unique_values
    
    if total_unique_combinations <= NUM_SWARMS:
        # CASE 1: Fewer unique combinations than swarms - use combinatorial approach
        var swarms_per_combo = NUM_SWARMS // total_unique_combinations
        var remainder = NUM_SWARMS % total_unique_combinations
        
        print("  Strategy: Distributing", swarms_per_combo, "swarms per combination +", remainder, "extra")
        
        # Generate all possible values for combinatorial approach
        var queen_values = List[Int]()
        var transition_values = List[Int]()
        var forager_values = List[Int]()
        var hypersensitivity_values = List[Int]()
        
        for i in range(queen_min, queen_max + 1):
            queen_values.append(i)
        for i in range(transition_min, transition_max + 1):
            transition_values.append(i)
        for i in range(forager_min, forager_max + 1):
            forager_values.append(i)
        for i in range(hypersensitivity_min, hypersensitivity_max + 1):
            hypersensitivity_values.append(i)
        
        # If no variation, ensure at least one value
        if len(queen_values) == 0:
            queen_values.append(queen_min)
        if len(transition_values) == 0:
            transition_values.append(transition_min)
        if len(forager_values) == 0:
            forager_values.append(forager_min)
        if len(hypersensitivity_values) == 0:
            hypersensitivity_values.append(hypersensitivity_min)
        
        var swarm_count = 0
        for queen_val in queen_values:
            for transition_val in transition_values:
                for forager_val in forager_values:
                    for hypersensitivity_val in hypersensitivity_values:
                        var num_for_this_combo = swarms_per_combo
                        if remainder > 0:
                            num_for_this_combo += 1
                            remainder -= 1
                        
                        for _ in range(num_for_this_combo):
                            if swarm_count >= NUM_SWARMS:
                                break
                            
                            var chromosome = Chromosome()
                            chromosome.allele1.append(queen_val)
                            chromosome.allele2.append(queen_val)
                            chromosome.allele1.append(transition_val)
                            chromosome.allele2.append(transition_val)
                            chromosome.allele1.append(forager_val)
                            chromosome.allele2.append(forager_val)
                            chromosome.allele1.append(hypersensitivity_val)
                            chromosome.allele2.append(hypersensitivity_val)
                            
                            swarms.append(generate_swarm_with_chromosome(chromosome))
                            swarm_count += 1
        
        # Fill any remaining slots with random values
        while swarm_count < NUM_SWARMS:
            var chromosome = Chromosome()
            chromosome.generate_alleles(queen_min, queen_max, transition_min, transition_max, forager_min, forager_max, hypersensitivity_min, hypersensitivity_max)
            swarms.append(generate_swarm_with_chromosome(chromosome))
            swarm_count += 1
            
    else:
        # CASE 2: More unique combinations than swarms - NEW APPROACH: equally distributed sampling
        print("  Strategy: Sampling", NUM_SWARMS, "equally distributed points from each allele range")
        
        # Generate equally distributed values for each allele
        var sampled_queen_values = generate_equally_distributed_values(queen_min, queen_max, NUM_SWARMS)
        var sampled_transition_values = generate_equally_distributed_values(transition_min, transition_max, NUM_SWARMS)
        var sampled_forager_values = generate_equally_distributed_values(forager_min, forager_max, NUM_SWARMS)
        var sampled_hypersensitivity_values = generate_equally_distributed_values(hypersensitivity_min, hypersensitivity_max, NUM_SWARMS)

        # Apply shuggling so the intersections are random
        shuffle(sampled_queen_values)
        shuffle(sampled_transition_values)
        shuffle(sampled_forager_values)
        shuffle(sampled_hypersensitivity_values)
        
        print("  Sampled queen values:", len(sampled_queen_values), "unique values")
        print("  Sampled nurse values:", len(sampled_transition_values), "unique values") 
        print("  Sampled forager values:", len(sampled_forager_values), "unique values")
        print("  Sampled hypersensitivity values:", len(sampled_hypersensitivity_values), "unique values")
        
        # Create swarms using the sampled values
        for swarm_idx in range(NUM_SWARMS):
            var chromosome = Chromosome()
            
            # Use the pre-sampled equally distributed values
            var queen_val = sampled_queen_values[swarm_idx]
            var transition_val = sampled_transition_values[swarm_idx]
            var forager_val = sampled_forager_values[swarm_idx]
            var hypersensitivity_val = sampled_hypersensitivity_values[swarm_idx]
            
            chromosome.allele1.append(queen_val)
            chromosome.allele2.append(queen_val)
            chromosome.allele1.append(transition_val)
            chromosome.allele2.append(transition_val)
            chromosome.allele1.append(forager_val)
            chromosome.allele2.append(forager_val)
            chromosome.allele1.append(hypersensitivity_val)
            chromosome.allele2.append(hypersensitivity_val)
            
            swarms.append(generate_swarm_with_chromosome(chromosome))
    
    print("Created", len(swarms), "swarms with optimized allele distribution")
    return swarms.copy()

fn reproduce_dead_swarm(mut dead_swarm: Swarm,
 parent1_chromosome: Chromosome,
 parent2_chromosome: Chromosome,
 mut rng: Xoshiro256PlusPlus):
    """Reproduce a dead swarm using genetic crossover from two parents."""
    
    # Create offspring chromosome through proper crossover
    var offspring_chromosome = Chromosome()
    
    # For each gene, randomly inherit one allele from each parent
    for gene_idx in range(NUM_GENES):
        # Inherit one allele from parent 1 (randomly choose allele1 or allele2)
        if rng.next_float64() < 0.5:
            offspring_chromosome.allele1.append(parent1_chromosome.allele1[gene_idx])
        else:
            offspring_chromosome.allele1.append(parent1_chromosome.allele2[gene_idx])
        
        # Inherit one allele from parent 2 (randomly choose allele1 or allele2)
        if rng.next_float64() < 0.5:
            offspring_chromosome.allele2.append(parent2_chromosome.allele1[gene_idx])
        else:
            offspring_chromosome.allele2.append(parent2_chromosome.allele2[gene_idx])

        
    # mutate on 5% of the value
    mutate_chromosome(offspring_chromosome, rng, mutation_chance=MUTATION_RATE)
    
    # Generate a new swarm with the mutated chromosome
    var new_swarm = generate_swarm_with_chromosome(offspring_chromosome)
    
    # Copy the new swarm's properties to the dead swarm
    dead_swarm.chromosome = new_swarm.chromosome.copy()
    dead_swarm.bees = new_swarm.bees.copy()
    dead_swarm.newborns_today = 0
    
    # Reset egg tracking
    for i in range(len(dead_swarm.eggs_by_day)):
        dead_swarm.eggs_by_day[i] = 0

fn force_trait_to_dominance(mut chromosome: Chromosome, gene_number: Int):
    """Choose allele which would have dominance and recessive mechanics (smaller value appear only if two alleles are recessive).
    Used for hypersensitivity."""
    var allele1 = chromosome.allele1[gene_number]
    var allele2 = chromosome.allele2[gene_number]

    if allele1 > allele2:
        chromosome.allele2[gene_number] = allele1
    elif allele2 > allele1:
        chromosome.allele1[gene_number] = allele2
    elif allele1 == allele2:
        pass


fn mutate_chromosome(mut chromosome: Chromosome,
    mut rng: Xoshiro256PlusPlus,
    mutation_chance: Float64 = MUTATION_RATE,
    min_mutation_amount: Int = 1,
    max_mutation_amount: Int = 1000,
    ):
    """Apply random mutation to chromosome with mutation_chance chance."""
    # Check if mutation occurs (mutation_chance)
    if rng.next_float64() < mutation_chance:
        print("Mutation triggered")
        
        # Safety check: ensure chromosome has proper structure
        if len(chromosome.allele1) != NUM_GENES or len(chromosome.allele2) != NUM_GENES:
            print("ERROR: Invalid chromosome structure in mutation")
            return
        
        # Randomly select which gene to mutate (0, 1, or 2) exclude hypersensitivity gene from mutation
        var gene_to_mutate = Int(random_si64(0, NUM_GENES-2))
        
        # Randomly select which allele to mutate
        var allele_to_mutate = Int(random_si64(0, 1))  # 0 or 1
        
        # Generate mutation value with log10 distribution (smaller values more likely)
        var log_min = log10(Float64(min_mutation_amount))
        var log_max = log10(Float64(max_mutation_amount))
        var u = rng.next_float64()  # Uniform [0, 1)
        var log_val = log_min + (log_max - log_min) * u
        var raw_mutation_value = pow(10.0, log_val)
        var mutation_value = Int(round(raw_mutation_value))
        
        # Clamp to bounds (in case of rounding edge cases)
        mutation_value = max(min_mutation_amount, min(max_mutation_amount, mutation_value))
        
        # Apply mutation
        var new_value = mutation_value
        
        # Apply the mutation safely
        if allele_to_mutate == 0:
            chromosome.allele1[gene_to_mutate] = new_value
        else:
            chromosome.allele2[gene_to_mutate] = new_value


fn check_swarm_death(mut swarm: Swarm,
 swarm_death_rate: Float64 = RANDOM_SWARM_DEATH):
    """Check if swarm dies randomly and clear all bees if it does."""
    if random_float64() < swarm_death_rate:
        # Swarm dies - clear all bees and reset eggs
        swarm.bees = List[Bee]()
        swarm.newborns_today = 0
        # Clear all eggs as well
        for i in range(len(swarm.eggs_by_day)):
            swarm.eggs_by_day[i] = 0

fn process_eggs_and_hatching(mut swarm: Swarm, pathogens: List[Pathogen],
 instant_nurses: Bool = True):
    """Handle daily egg laying and hatching.
    Args:
        instant_nurses: If True, nurses are added immediately without n-day delay.
    """
    # Skip processing if swarm is dead (no bees)
    if len(swarm.bees) == 0:
        swarm.newborns_today = 0
        return
    
    # Reset today's newborns count
    swarm.newborns_today = 0
    
    # Check if queen exists and is fertile
    var has_fertile_queen = is_queen_alive_and_fertile(swarm, pathogens)
    var current_bee_count = len(swarm.bees)
    var space_available = MAX_BEES_PER_SWARM - current_bee_count
    var current_fitness = calculate_fitness(swarm, pathogens)
    
    if instant_nurses:
        # INSTANT NURSES MODE: Add nurses immediately if there's space AND queen exists
        if space_available > 0 and has_fertile_queen:
            var nurses_to_add = min(EGGS_PER_DAY, space_available)

            # check does fitness allow to add all of the nurses
            if current_fitness >= nurses_to_add * 1000:
                # Add new nurse bees immediately
                for _ in range(nurses_to_add):
                    swarm.bees.append(Bee(caste_type = NURSE, base_longevity = generate_base_longevity(NURSE, swarm.chromosome)))
                    swarm.newborns_today += 1
            else:
                nurses_to_add = current_fitness // 1000
                # Add new nurse bees immediately
                for _ in range(nurses_to_add):
                    swarm.bees.append(Bee(caste_type = NURSE, base_longevity = generate_base_longevity(NURSE, swarm.chromosome)))
                    swarm.newborns_today += 1

            # print("Nurses added: ", nurses_to_add)
    else:
        # ORIGINAL MODE: n-day development delay
        # Hatch eggs that are ready (from n days ago) only if there's space
        var eggs_to_hatch = swarm.eggs_by_day[EGG_DEVELOPMENT_DAYS - 1]
        
        if space_available > 0:
            var nurses_to_add = min(eggs_to_hatch, space_available)
            
            # Add new nurse bees and count them as newborns
            for _ in range(nurses_to_add):
                swarm.bees.append(Bee(caste_type = NURSE, base_longevity = generate_base_longevity(NURSE, swarm.chromosome)))
                swarm.newborns_today += 1
        
        # Shift egg development forward by one day
        for i in range(EGG_DEVELOPMENT_DAYS - 1, 0, -1):
            swarm.eggs_by_day[i] = swarm.eggs_by_day[i - 1]
        # Clear the first day after shifting
        swarm.eggs_by_day[0] = 0
    
    # Lay new eggs (only if there's a fertile queen and space for future growth)
    var should_lay_eggs = (has_fertile_queen and len(swarm.bees) < MAX_BEES_PER_SWARM)
    
    if should_lay_eggs:
        if instant_nurses:
            # In instant mode, we've already added nurses, so lay eggs for future
            # But only if we have space after the nurses we just added
            if space_available >= EGGS_PER_DAY:
                swarm.eggs_by_day[0] = EGGS_PER_DAY
            else:
                swarm.eggs_by_day[0] = 0
        else:
            # In normal mode, eggs go into the development pipeline
            swarm.eggs_by_day[0] = EGGS_PER_DAY
    else:
        swarm.eggs_by_day[0] = 0

fn age_and_kill_bees(mut swarm: Swarm, pathogens: List[Pathogen],
 mut simulation_data: SimulationData,
 mut rng: Xoshiro256PlusPlus
) -> Tuple[Int, Int, Int, Int]:
    """Age all bees by 1 day, handle transitions, remove those that die.
    Apply hypersensitivity death on infected nurses and foragers.
     Also add swarm dead data to simulation_data."""
    # Skip if swarm is dead
    if len(swarm.bees) == 0:
        return (0, 0, 0, 0)
        
    var survivors = List[Bee]()

    var dead_by_aging = 0
    var dead_by_non_aging = 0
    var dead_by_pathogen = 0
    var dead_by_hypersensitivity = 0

    var hypersensitivity_probability: Float64 = Float64(swarm.chromosome.get_trait_value(3)) / 100 # used to check probability of an infected worker to die from hypersensitivity
    var transition_age: Int = Int(swarm.chromosome.get_trait_value(TRANSITION_AGE_GENE))

    for i in range(len(swarm.bees)):
        var ref bee = swarm.bees[i]
        # Age the bee
        bee.age += 1
        
        # Check for nurse to forager transition
        if bee.caste_type == NURSE and bee.age >= transition_age:
            # Transition nurse to forager
            bee.caste_type = FORAGER
            # Reset age
            bee.age = 0 
            # Update longevity for new forager role
            bee.base_longevity = generate_base_longevity(FORAGER, swarm.chromosome)
        
        # Check death causes in order of priority
        var bee_died = False
        
        # 1. Check if bee dies from old age
        if bee.age >= bee.base_longevity:
            bee_died = True
            dead_by_aging += 1 # Track aging deaths by caste, nurses do not die from aging
        
        # 2. Check if forager dies from random death (only if not already dead)
        if not bee_died and bee.caste_type == FORAGER:
            if rng.next_float() < FORAGER_DEATH_RATE:
                bee_died = True
                dead_by_non_aging += 1
                
        # 3. Check if bee dies from pathogen-induced death (only if not already dead)
        if not bee_died and bee.infected != 0:  # Check if bee has any infection
            # Check each pathogen that the bee is infected with
            for pathogen_idx in range(len(pathogens)):
                if pathogens[pathogen_idx].increases_death_chance:
                    var pathogen_bit = UInt8(1) << pathogen_idx
                    var pathogen = pathogens[pathogen_idx].copy()
                        
                    if (bee.infected & pathogen_bit) != 0: # Check if bee is infected by this specific pathogen
                        if rng.next_float() < pathogen.death_probability:
                            bee_died = True
                            dead_by_pathogen += 1
                            break  # Bee died, no need to check other pathogens

        # 4. Check if infected worker dies by hypersensitivity
        if not bee_died and bee.infected != 0 and (bee.caste_type == NURSE or bee.caste_type == FORAGER):
            if hypersensitivity_probability > rng.next_float():
                bee_died = True
                dead_by_hypersensitivity += 1

        # 5 Check if queen is sterile and is the only bee in a hive, kill her
        if bee.caste_type == QUEEN and len(swarm.bees) == 1 and not is_queen_alive_and_fertile(swarm, pathogens):
            bee_died = True
            dead_by_non_aging += 1 # Count sterile queen death as non-aging death

        # Keep bee if it doesn't die from any cause
        if not bee_died:
            survivors.append(bee.copy())
    
    # Replace original swarm with survivors
    swarm.bees = survivors.copy()

    return (dead_by_aging, dead_by_non_aging, dead_by_pathogen, dead_by_hypersensitivity)


fn attempt_cure(mut swarm: Swarm, pathogens: List[Pathogen],
 mut rng: Xoshiro256PlusPlus):
    """Parallel bee treatment - ACTUALLY WORKS with bitmask logic."""
    var num_bees = len(swarm.bees)
    
    @parameter
    fn treat_bee(i: Int):
        # This works because parallelize handles closures differently
        ref bee = swarm.bees[i]
        
        # Check each pathogen for potential cure/disappearance
        for pathogen_idx in range(len(pathogens)):
            var pathogen_bit = UInt8(1) << pathogen_idx
            
            # Check if bee is infected by this pathogen
            if (bee.infected & pathogen_bit) != 0:
                # Check if pathogen disappears with some probability
                if rng.next_float() < pathogens[pathogen_idx].disappearing_probability:
                    # Clear infection for this pathogen
                    bee.infected &= ~pathogen_bit
                    
                    # Set immunity for this pathogen (if applicable)
                    # Uncomment if pathogens grant immunity after disappearance:
                    bee.immune |= pathogen_bit
                    
                    # Note: With bitmask, we don't need to break since bee can have multiple infections
                    # If you want only one cure per bee per day, uncomment:
                    # break
    
    # This WILL work with parallelize!
    parallelize[treat_bee](num_bees)

