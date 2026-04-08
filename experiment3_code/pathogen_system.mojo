from types import *
from math import floor
from algorithm import parallelize
from random import randint, rand
from randomizer import Xoshiro256PlusPlus

comptime OUT_OF_SWARM_FORAGER_INFECTION_PROB: Float64 = 0.0001

fn count_working_workers(swarm: Swarm, pathogens: List[Pathogen]) -> Tuple[Int, Int]:
    """Count working nurses and foragers that are able to work (not disabled by pathogens). Returns tuple of two numbers."""
    var working_nurses = 0
    var working_foragers = 0

    for bee in swarm.bees:
        if bee.caste_type == NURSE:
            var can_work = True
            
            # Check if bee is infected with any pathogen that disables work
            if bee.infected != 0:  # Bee has at least one infection
                # Check each pathogen bit
                for pathogen_idx in range(len(pathogens)):
                    # Check if bee is infected by this pathogen
                    if (bee.infected & (1 << pathogen_idx)) != 0:
                        var pathogen = pathogens[pathogen_idx].copy()
                        if pathogen.disables_worker_work:
                            can_work = False
                            break  # No need to check other pathogens
            
            if can_work:
                working_nurses += 1

        elif bee.caste_type == FORAGER:
            var can_work = True
            
            # Check if bee is infected with any pathogen that disables work
            if bee.infected != 0:  # Bee has at least one infection
                # Check each pathogen bit
                for pathogen_idx in range(len(pathogens)):
                    # Check if bee is infected by this pathogen
                    if (bee.infected & (1 << pathogen_idx)) != 0:
                        var pathogen = pathogens[pathogen_idx].copy()
                        if pathogen.disables_worker_work:
                            can_work = False
                            break  # No need to check other pathogens
            
            if can_work:
                working_foragers += 1
    
    return Tuple(working_nurses, working_foragers)

fn is_queen_alive_and_fertile(swarm: Swarm, pathogens: List[Pathogen]) -> Bool:
    """Check if there's a living, fertile queen."""
    for bee in swarm.bees:
        if bee.caste_type == QUEEN:
            # Queen is alive, now check if fertile
            
            # Check if queen is infected with any pathogen that causes infertility
            if bee.infected != 0:  # Queen has at least one infection
                # Check each pathogen bit
                for pathogen_idx in range(len(pathogens)):
                    # Check if queen is infected by this pathogen
                    if (bee.infected & (1 << pathogen_idx)) != 0:
                        var pathogen = pathogens[pathogen_idx].copy()
                        if pathogen.makes_queen_infertile:
                            return False  # Queen is infertile
            
            return True  # Queen is alive and fertile
    return False  # No queen found

fn is_queen_fertile(swarm: Swarm, pathogens: List[Pathogen]) -> Bool:
    """Check if the queen can lay eggs (not infected with infertility-causing pathogen)."""
    for bee in swarm.bees:
        if bee.caste_type == QUEEN and bee.infected != 0:
            for pathogen_idx in range(len(pathogens)):
                var pathogen = pathogens[pathogen_idx].copy()
                if pathogen.makes_queen_infertile:
                    return False
    return True

fn infect_foragers(mut swarm: Swarm, pathogens: List[Pathogen],
 mut rng: Xoshiro256PlusPlus):
    """Randomly infect foragers with pathogens, respecting immunity."""
    # Skip if swarm is dead
    if len(swarm.bees) == 0:
        return
        
    for i in range(len(swarm.bees)):
        var ref bee = swarm.bees[i]
        # Only foragers can be initially infected
        if bee.caste_type == FORAGER:
            # Check each pathogen for potential infection
            for pathogen_idx in range(len(pathogens)):

                if pathogens[pathogen_idx].infectivity > 0.0:
                    # Check immunity first
                    var pathogen_bit = UInt8(1) << pathogen_idx
                    var is_immune = (bee.immune & pathogen_bit) != 0
                    
                    if not is_immune:
                        # Check if already infected
                        var already_infected = (bee.infected & pathogen_bit) != 0
                        
                        if not already_infected:
                            if rng.next_float() < OUT_OF_SWARM_FORAGER_INFECTION_PROB:  # 0.01% chance per day per pathogen
                                    # Set the infection bit for this pathogen
                                    bee.infected |= pathogen_bit
       

fn spread_infection(mut swarm: Swarm, pathogens: List[Pathogen],
 mut rng: Xoshiro256PlusPlus,
 ):
    """Spread infection based on number of infected bees and beta parameter. Each pathogen is processed one by one."""
    # Skip if swarm is dead
    var num_bees = len(swarm.bees)
    if num_bees == 0:
        return
    
    # Initialize count list for all pathogens
    var infected_counts = List[Int]()
    for _ in range(len(pathogens)):
        infected_counts.append(0)
    
    # Single pass counting of infected bees for all pathogens
    for bee in swarm.bees:
        if bee.infected != 0:
            for pathogen_idx in range(len(pathogens)):
                var pathogen_bit = UInt8(1) << pathogen_idx
                if (bee.infected & pathogen_bit) != 0:
                    infected_counts[pathogen_idx] += 1
    
    # Loop over each pathogen to apply infections
    for pathogen_idx in range(len(pathogens)):
        var number_of_infected_bees_by_this_pathogen = infected_counts[pathogen_idx]
        
        if number_of_infected_bees_by_this_pathogen == 0:
            continue
        
        # Calculate number of infection shots
        var infection_shots = calculate_infection_shots(
            number_of_infected_bees_by_n_pathogen=number_of_infected_bees_by_this_pathogen,
            infectivity=pathogens[pathogen_idx].infectivity)
        
        var pathogen_bit = UInt8(1) << pathogen_idx

        # Perform infection shots with original random logic
        for _ in range(infection_shots):
            var target_index = rng.next_range(0, MAX_BEES_PER_SWARM)
            
            if target_index < num_bees:
                var ref bee = swarm.bees[target_index]
                
                # Check if bee is already immune to this pathogen
                var is_immune = (bee.immune & pathogen_bit) != 0
                
                if not is_immune:
                    # Set infection bit for this pathogen
                    bee.infected |= pathogen_bit
                    # Note: With bitmask, bee can be infected by multiple pathogens

              
fn calculate_infection_shots(number_of_infected_bees_by_n_pathogen: Int, infectivity: Float64) -> Int:
    """Caclulate infection shots based on number of already infected and pathogen infectivity,  used in spread_infection function."""
    var shot_number = number_of_infected_bees_by_n_pathogen * infectivity
    var shot_number_residue = shot_number - floor(shot_number)

    if  random_float64() > shot_number_residue:    #  random_float64()
        result =  Int(floor(shot_number))
    else:
        result =  Int(floor(shot_number) + 1)

    return result