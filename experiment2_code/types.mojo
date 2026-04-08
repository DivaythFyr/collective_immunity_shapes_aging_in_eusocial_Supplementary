from memory import memset_zero
from random import random_float64, random_si64
from analysis import calculate_fitness
from algorithm import parallelize
from main import *

# Bee caste enumeration (removed DRONE)
comptime QUEEN: Int = 0
comptime NURSE: Int = 1
comptime FORAGER: Int = 2

# Constants
comptime MAX_BEES_PER_SWARM: Int = 2001
comptime NUM_SWARMS: Int = 300
comptime EGGS_PER_DAY: Int = Int((MAX_BEES_PER_SWARM - 1) // 100 * 10)
comptime EGG_DEVELOPMENT_DAYS  = 1

# Initial distribution per swarm (without drones, redistributed)
comptime INITIAL_QUEENS: Int = 1
comptime INITIAL_NURSES: Int = (MAX_BEES_PER_SWARM-1) // 100 * 80
comptime INITIAL_FORAGERS: Int = (MAX_BEES_PER_SWARM-1) // 100 * 20

# How much contribution provide worker classes to fitness calculation
comptime MAX_NURSE_CONTRIBUTION: Int = (MAX_BEES_PER_SWARM-1) // 100 * 80
comptime MAX_FORAGER_CONTRIBUTION:  Int = (MAX_BEES_PER_SWARM-1) // 100 * 20
comptime MAX_FITNESS: Int = MAX_NURSE_CONTRIBUTION * MAX_FORAGER_CONTRIBUTION

# Gene indices for chromosome
comptime QUEEN_AGE_GENE: Int = 0
comptime TRANSITION_AGE_GENE: Int = 1
comptime FORAGER_AGE_GENE: Int = 2
comptime HYPERSENSITIVITY_GENE: Int = 3
comptime NUM_GENES: Int = 4

struct Pathogen(Copyable, Movable):
    """Represents a pathogen that can infect bees."""
    var ID: Int
    var infectivity: Float64  # Probability of infection per contact (0.0 to 1.0)
    var makes_queen_infertile: Bool  # Does it make queens infertile when infected
    var disables_worker_work: Bool # Does it make nurses unable to work when infected
    var increases_death_chance: Bool  # Does it increase chance of death
    var death_probability: Float64  # Daily death probability if increases_death_chance is True (0.1 to 0.9)
    var disappearing_probability: Float64 # chance of a bee to be cured by this pathogen
    
    fn __init__(out self, ID: Int, infectivity: Float64, makes_queen_infertile: Bool, disables_worker_work: Bool, increases_death_chance: Bool, death_probability: Float64, disappearing_probability: Float64):
        self.ID = ID
        self.infectivity = infectivity
        self.makes_queen_infertile = makes_queen_infertile
        self.disables_worker_work = disables_worker_work
        self.increases_death_chance = increases_death_chance
        self.death_probability = death_probability
        self.disappearing_probability = disappearing_probability

struct Chromosome(Copyable, Movable):
    """Represents genetic information with two alleles for each gene."""
    var allele1: List[Int]  # First allele for each gene
    var allele2: List[Int]  # Second allele for each gene
    
    fn __init__(out self):
        self.allele1 = List[Int]()
        self.allele2 = List[Int]()

    fn generate_alleles(mut self, queen_min: Int, queen_max: Int, transition_min: Int, transition_max: Int, forager_min: Int, forager_max: Int, hypersensitivity_min: Int, hypersensitivity_max: Int):
        """Fill a chromosome with specific allele ranges. Two alleles are lists that are filled with a random value from defined min to max."""
        # Queen base longevity gene: use specified range
        self.allele1.append(queen_min + Int(random_si64(0, queen_max - queen_min)))
        self.allele2.append(queen_min + Int(random_si64(0, queen_max - queen_min)))
        
        # Transition base longevity gene: use specified range
        self.allele1.append(transition_min + Int(random_si64(0, transition_max - transition_min)))
        self.allele2.append(transition_min + Int(random_si64(0, transition_max - transition_min)))
        
        # Forager base longevity gene: use specified range
        self.allele1.append(forager_min + Int(random_si64(0, forager_max - forager_min)))
        self.allele2.append(forager_min + Int(random_si64(0, forager_max - forager_min)))

        # Hypersensitivity gene: use specified range
        self.allele1.append(hypersensitivity_min + Int(random_si64(0, hypersensitivity_max - hypersensitivity_min)))
        self.allele2.append(hypersensitivity_min + Int(random_si64(0, hypersensitivity_max - hypersensitivity_min)))
    
    fn get_trait_value(self, gene_index: Int) -> Int:
        """Get the actual trait value as average of two alleles."""
        if gene_index == 3:
            return max(self.allele1[gene_index], self.allele2[gene_index])
        else:
            return (self.allele1[gene_index] + self.allele2[gene_index]) // 2

struct Bee(Copyable, Movable):
    """Represents a bee caste with its properties."""
    var caste_type: Int  # 0=Queen, 1=Nurse, 2=Forager
    var base_longevity: Int  # Base longevity for this caste
    var age: Int  # Current age in days
    var infected: UInt8
    var immune: UInt8
    
    fn __init__(out self, caste_type: Int, base_longevity: Int, age: Int = 0):
        self.caste_type = caste_type
        self.base_longevity = base_longevity
        self.age = age
        self.infected = UInt8(0)
        self.immune = UInt8(0)

struct Swarm(Copyable, Movable):
    """Represents a swarm with its bees and swarm-specific parameters."""
    var bees: List[Bee]
    var chromosome: Chromosome  # Genetic information for the swarm
    var eggs_by_day: List[Int]  # Eggs that will hatch in 1, 2, 3, 4, 5, 6 days
    var newborns_today: Int  # Number of bees that hatched today
    
    fn __init__(out self):
        self.bees = List[Bee]()
        self.chromosome = Chromosome()
        self.eggs_by_day = List[Int]()
        self.newborns_today = 0
        
        # Initialize egg tracking (n days until hatching)
        for _ in range(EGG_DEVELOPMENT_DAYS):
            self.eggs_by_day.append(0)

struct SwarmFitness(Copyable, Movable):
    """Stores swarm index and its fitness for sorting."""
    var swarm_index: Int
    var fitness: Int
    
    fn __init__(out self, swarm_index: Int, fitness: Int):
        self.swarm_index = swarm_index
        self.fitness = fitness

struct DailySwarmData(Copyable, Movable):
    """Stores daily population data for one swarm."""
    var queens: Int
    var nurses: Int
    var foragers: Int
    
    fn __init__(out self, queens: Int, nurses: Int, foragers: Int):
        self.queens = queens
        self.nurses = nurses
        self.foragers = foragers

struct DailyAlleleData(Copyable, Movable):
    """Stores daily allele data for one swarm."""
    var swarm_id: Int
    var queen_age_allele: Int
    var transition_age_allele: Int
    var forager_age_allele: Int
    var hypersensitivity_allele: Int
    var fitness: Int
    
    fn __init__(out self, swarm_id: Int, queen_age_allele: Int, transition_age_allele: Int, forager_age_allele: Int, hypersensitivity_allele: Int, fitness: Int):
        self.swarm_id = swarm_id
        self.queen_age_allele = queen_age_allele
        self.transition_age_allele = transition_age_allele
        self.forager_age_allele = forager_age_allele
        self.hypersensitivity_allele = hypersensitivity_allele
        self.fitness = fitness

struct SimulationData(Copyable, Movable):
    """Stores all simulation data for visualization."""
    var daily_data: List[List[DailySwarmData]]
    var daily_allele_data: List[List[DailyAlleleData]]
    var daily_infection_data: List[List[DailyInfectionData]]
    var days: Int
    var run_id: String
    # Death tracking fields:
    
    # ← ADD SWARM DEATH TRACKING:
    var total_swarms_died: Int           # Total number of swarms that have died during simulation

    # Keep in mind that these parameters are changed in age_and_kill_bees function, not from add_day_data method
    # Each element in the outer list is day. Each element in the inner list is dead quantity for one swarm
    var daily_aging_dead: List[List[Int]]
    var daily_non_aging_dead: List[List[Int]]
    var daily_pathogen_dead: List[List[Int]]
    var daily_hypersensitivity_dead: List[List[Int]]

    var infected_by_pathogen: List[List[Int]]
    var immune_to_pathogen: List[List[Int]]
    
    fn __init__(out self, run_id: String = "run"):
        self.daily_data = List[List[DailySwarmData]]()
        self.daily_allele_data = List[List[DailyAlleleData]]()
        self.daily_infection_data = List[List[DailyInfectionData]]()
        self.days = 0
        self.run_id = run_id
        
        # ← INITIALIZE SWARM DEATH COUNTER:
        self.total_swarms_died = 0

        self.daily_aging_dead = List[List[Int]]()
        self.daily_non_aging_dead = List[List[Int]]()
        self.daily_pathogen_dead = List[List[Int]]()
        self.daily_hypersensitivity_dead = List[List[Int]]()

        self.infected_by_pathogen = List[List[Int]]()
        self.immune_to_pathogen = List[List[Int]]()

    fn add_day_data(mut self, swarms: List[Swarm], pathogens: List[Pathogen]):
        """Add data for one day (parallelized per-swarm processing)."""
        var day_data = List[DailySwarmData](capacity=len(swarms))
        var day_allele_data = List[DailyAlleleData](capacity=len(swarms))
        var day_infection_data = List[DailyInfectionData](capacity=len(swarms))

        # Initialize per-swarm result lists
        for _ in range(len(swarms)):
            day_data.append(DailySwarmData(0, 0, 0))
            day_allele_data.append(DailyAlleleData(0, 0, 0, 0, 0, 0))
            day_infection_data.append(DailyInfectionData(0, 0, 0))

        # Initialize infection counters for all pathogens (thread-safe aggregation after parallel)
        var infected_by_pathogen_day = List[Int]()
        var immune_to_pathogen_day = List[Int]()
        for _ in range(len(pathogens)):
            infected_by_pathogen_day.append(0)
            immune_to_pathogen_day.append(0)

        # Parallel per-swarm processing
        @parameter
        fn process_swarm(i: Int):
            var swarm = swarms[i].copy()
            var counts = count_bees_by_caste(swarm)
            day_data[i] = DailySwarmData(counts[0], counts[1], counts[2])

            swarm.newborns_today = 0
            var fitness = calculate_fitness(swarm, pathogens)

            var queen_age_allele = swarm.chromosome.get_trait_value(QUEEN_AGE_GENE)
            var transition_age_allele = swarm.chromosome.get_trait_value(TRANSITION_AGE_GENE)
            var forager_age_allele = swarm.chromosome.get_trait_value(FORAGER_AGE_GENE)
            var hypersensitivity_allele = swarm.chromosome.get_trait_value(HYPERSENSITIVITY_GENE)

            day_allele_data[i] = DailyAlleleData(i, queen_age_allele, transition_age_allele, forager_age_allele, hypersensitivity_allele, fitness)

            var infection_counts = count_infected_bees_by_caste(swarm)
            day_infection_data[i] = DailyInfectionData(infection_counts[0], infection_counts[1], infection_counts[2])

        parallelize[process_swarm](len(swarms))

        # Serial pathogen infection/immunity counting (must be done after parallel)
        for i in range(len(swarms)):
            var swarm = swarms[i].copy()
            for bee in swarm.bees:
                for pathogen_idx in range(len(pathogens)):
                    var pathogen_bit = UInt8(1) << pathogen_idx
                    if (bee.infected & pathogen_bit) != 0:
                        infected_by_pathogen_day[pathogen_idx] += 1
                    if (bee.immune & pathogen_bit) != 0:
                        immune_to_pathogen_day[pathogen_idx] += 1

        self.daily_data.append(day_data.copy())
        self.daily_allele_data.append(day_allele_data.copy())
        self.daily_infection_data.append(day_infection_data.copy())
        self.days += 1

        self.infected_by_pathogen.append(infected_by_pathogen_day^)
        self.immune_to_pathogen.append(immune_to_pathogen_day^)


    fn get_pathogen_infection_immunity_pairs(self, day: Int, pathogens: List[Pathogen]) -> List[Tuple[Int, Int]]:
        """Get pairs of (infected_count, immune_count) for all pathogens for a specific day.
        Args:
            day: The day index to retrieve data for.
        Returns:
            List of tuples where each tuple is (infected_count, immune_count) for a pathogen
            The list length equals the number of pathogens.
        """
        var result = List[Tuple[Int, Int]]()
        
        # Check if day is valid
        if day < 0 or day >= self.days:
            print("Warning: Invalid day", day, "in get_pathogen_infection_immunity_pairs")
            return result^
        
        # Check if we have data for this day
        if day >= len(self.infected_by_pathogen) or day >= len(self.immune_to_pathogen):
            print("Warning: No pathogen data for day", day)
            return result^
        
        var infected_day = self.infected_by_pathogen[day].copy()
        var immune_day = self.immune_to_pathogen[day].copy()
        
        # Number of pathogens might differ between lists, take the smaller size
        var num_pathogens = len(pathogens)
        
        for pathogen_idx in range(num_pathogens):
            var infected_count = Int(infected_day[pathogen_idx])
            var immune_count = Int(immune_day[pathogen_idx])
            result.append((infected_count, immune_count))
        
        return result^

fn count_bees_by_caste(swarm: Swarm) -> Tuple[Int, Int, Int]:
    """Count bees by caste type. Returns (queens, nurses, foragers)."""
    var queens = 0
    var nurses = 0
    var foragers = 0
    
    for bee in swarm.bees:
        if bee.caste_type == QUEEN:
            queens += 1
        elif bee.caste_type == NURSE:
            nurses += 1
        elif bee.caste_type == FORAGER:
            foragers += 1
    
    return (queens, nurses, foragers)

struct DailyInfectionData(Copyable, Movable):
    """Stores daily infection data for one swarm."""
    var infected_queens: Int
    var infected_nurses: Int
    var infected_foragers: Int
    
    fn __init__(out self, infected_queens: Int, infected_nurses: Int, infected_foragers: Int):
        self.infected_queens = infected_queens
        self.infected_nurses = infected_nurses
        self.infected_foragers = infected_foragers

fn count_infected_bees_by_caste(swarm: Swarm) -> Tuple[Int, Int, Int]:
    """Count infected bees by caste type. Returns (infected_queens, infected_nurses, infected_foragers)."""
    var infected_queens = 0
    var infected_nurses = 0
    var infected_foragers = 0
    
    for bee in swarm.bees:
        if bee.infected != 0:
            if bee.caste_type == QUEEN:
                infected_queens += 1
            elif bee.caste_type == NURSE:
                infected_nurses += 1
            elif bee.caste_type == FORAGER:
                infected_foragers += 1
    
    return (infected_queens, infected_nurses, infected_foragers)