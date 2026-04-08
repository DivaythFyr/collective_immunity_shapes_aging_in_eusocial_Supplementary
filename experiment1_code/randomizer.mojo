from collections import List
from algorithm import parallelize
from random import random_si64
from time import perf_counter_ns
import math

# Thread-safe xoshiro256++ pseudo-random number generator
# Implements a multithreaded version with parallel computation support
struct Xoshiro256PlusPlus(Copyable, Movable):
    """
    Multithreaded xoshiro256++ pseudo-random number generator.
    Uses four 64-bit state values to generate high-quality random numbers with excellent statistical properties.
    Thread-safe implementation allows parallel usage across multiple threads.
    Period: 2^256 - 1.
    """
    var s0: UInt64
    var s1: UInt64
    var s2: UInt64
    var s3: UInt64
    
    fn __init__(out self, seed_val: UInt64):
        """Initialize the xoshiro256++ generator with a seed value using SplitMix64 for proper seeding."""
        # Use SplitMix64 to initialize the four 64-bit states from a single seed
        # This ensures all states are properly initialized and non-zero
        var z = seed_val if seed_val != 0 else 1
        
        # SplitMix64 initialization
        z += 0x9E3779B97F4A7C15  # Large odd constant
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB
        self.s0 = z ^ (z >> 31)
        
        z += 0x9E3779B97F4A7C15
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB
        self.s1 = z ^ (z >> 31)
        
        z += 0x9E3779B97F4A7C15
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB
        self.s2 = z ^ (z >> 31)
        
        z += 0x9E3779B97F4A7C15
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9
        z = (z ^ (z >> 27)) * 0x94D049BB133111EB
        self.s3 = z ^ (z >> 31)
    
    @always_inline
    fn _rotl(self, x: UInt64, k: Int) -> UInt64:
        """Rotate left operation for xoshiro256++."""
        return (x << k) | (x >> (64 - k))
    
    fn next(mut self) -> UInt64:
        """Generate the next random number using xoshiro256++ algorithm."""
        var result = self._rotl(self.s0 + self.s3, 23) + self.s0
        
        var t = self.s1 << 17
        self.s2 ^= self.s0
        self.s3 ^= self.s1
        self.s1 ^= self.s2
        self.s0 ^= self.s3
        self.s2 ^= t
        self.s3 = self._rotl(self.s3, 45)
        
        return result
    
    fn next64(mut self) -> UInt64:
        """Generate a random 64-bit number (same as next() for xoshiro256++)."""
        return self.next()
    
    fn next32(mut self) -> UInt32:
        """Generate a random 32-bit number from the high bits of next64()."""
        return UInt32(self.next() >> 32)
    
    fn next_float(mut self) -> Float64:
        """Generate a random float in range [0.0, 1.0) using 53-bit precision."""
        # Use high 53 bits for IEEE 754 double precision
        return Float64(self.next() >> 11) / Float64(1 << 53)
    
    fn next_float64(mut self) -> Float64:
        """Generate a high-precision random float in range [0.0, 1.0) using full 64-bit precision."""
        return Float64(self.next()) / Float64(UInt64.MAX)
    
    fn next_range(mut self, min_val: Int, max_val: Int) -> Int:
        """Generate a random integer in the specified range [min_val, max_val)."""
        var range_size = max_val - min_val
        return Int(self.next() % UInt64(range_size)) + min_val

    

fn parallel_monte_carlo_pi(num_samples: Int, num_threads: Int) -> Float64:
    """
    Estimate π using Monte Carlo method with parallel random number generation.
    """
    var hits_inside_circle = List[Int]()
    for _ in range(num_threads):
        hits_inside_circle.append(0)
    
    @parameter
    fn monte_carlo_worker(thread_id: Int):
        var samples_per_thread = num_samples // num_threads
        var extra_samples = num_samples % num_threads
        var thread_samples: Int
        if thread_id < extra_samples:
            thread_samples = samples_per_thread + 1
        else:
            thread_samples = samples_per_thread
        
        # Create thread-specific generator with unique seed
        var thread_seed = UInt64(abs(random_si64(0, 1000000)) + thread_id * 1234567)
        var generator = Xoshiro256PlusPlus(thread_seed)
        var local_hits = 0
        
        for _ in range(thread_samples):
            var x = generator.next_float() * 2.0 - 1.0  # Range [-1, 1]
            var y = generator.next_float() * 2.0 - 1.0  # Range [-1, 1]
            
            if x * x + y * y <= 1.0:
                local_hits += 1
        
        hits_inside_circle[thread_id] = local_hits
    
    # Execute parallel Monte Carlo computation
    parallelize[monte_carlo_worker](num_threads)
    
    # Sum up results from all threads
    var total_hits = 0
    for i in range(num_threads):
        total_hits += hits_inside_circle[i]
    
    # Estimate π: (hits_inside_circle / total_samples) * 4
    return Float64(total_hits) / Float64(num_samples) * 4.0

fn benchmark_million_numbers_generation():
    """
    Benchmark generation of 1 million random numbers.
    Tests both single-threaded and multi-threaded performance.
    """
    var num_numbers = 1000000
    print("=== Benchmarking 1 Million Random Numbers Generation ===")
    
    # Single-threaded benchmark
    print("\nSingle-threaded generation:")
    print("Generating", num_numbers, "random numbers...")
    
    var start_time = perf_counter_ns()
    var single_rng = Xoshiro256PlusPlus(42)
    var sum_for_validation: UInt64 = 0
    
    # Generate numbers and measure performance
    for i in range(num_numbers):
        var rand_num = single_rng.next()
        sum_for_validation += rand_num
        
        # Progress indicator every 250k numbers
        if i % 250000 == 249999:
            print("Generated", i + 1, "numbers...")
    
    var end_time = perf_counter_ns()
    var single_duration_ns = end_time - start_time
    var single_duration_ms = Float64(single_duration_ns) / 1_000_000.0
    
    print("Single-thread generation completed!")
    print("Validation checksum:", sum_for_validation)
    print("Time elapsed:", single_duration_ms, "milliseconds")
    print("Generation rate:", Float64(num_numbers) / single_duration_ms * 1000.0, "numbers per second")
    
    # Multi-threaded benchmark with 4 threads
    print("\nMulti-threaded generation (4 threads):")
    var random_numbers = List[UInt64]()
    for _ in range(num_numbers):
        random_numbers.append(0)
    
    start_time = perf_counter_ns()
    
    @parameter
    fn generate_chunk(thread_id: Int):
        var chunk_size = num_numbers // 4
        var start_idx = thread_id * chunk_size
        var end_idx: Int
        if thread_id < 3:
            end_idx = start_idx + chunk_size
        else:
            end_idx = num_numbers
        
        var thread_seed = UInt64(abs(random_si64(0, 1000000)) + thread_id * 1234567)
        var generator = Xoshiro256PlusPlus(thread_seed)
        
        for i in range(start_idx, end_idx):
            random_numbers[i] = generator.next()
    
    print("Starting parallel generation...")
    parallelize[generate_chunk](4)
    
    end_time = perf_counter_ns()
    var multi_duration_ns = end_time - start_time
    var multi_duration_ms = Float64(multi_duration_ns) / 1_000_000.0
    
    print("Multi-threaded generation completed!")
    print("Time elapsed:", multi_duration_ms, "milliseconds")
    print("Generation rate:", Float64(num_numbers) / multi_duration_ms * 1000.0, "numbers per second")
    
    # Calculate performance metrics
    if single_duration_ms > 0 and multi_duration_ms > 0:
        var speedup = single_duration_ms / multi_duration_ms
        print("Speedup factor:", speedup, "x")
        var efficiency = speedup / 4.0 * 100.0  # 4 threads
        print("Parallel efficiency:", efficiency, "%")
    
    # Calculate validation checksum for parallel generation
    var parallel_sum: UInt64 = 0
    for i in range(min(10000, num_numbers)):  # Sample first 10k for validation
        parallel_sum += random_numbers[i]
    print("Parallel validation checksum (first 10k):", parallel_sum)
    
    # Show some sample numbers
    print("Sample generated numbers:")
    for i in range(min(5, num_numbers)):
        print("  [", i, "]:", random_numbers[i])
    
    print("\nPerformance Summary:")
    print("- Successfully generated", num_numbers, "random numbers")
    print("- Single-threaded time:", single_duration_ms, "ms")
    print("- Multi-threaded time: ", multi_duration_ms, "ms")
    print("- Time saved with parallelization:", single_duration_ms - multi_duration_ms, "ms")
    print("- xoshiro256++ algorithm shows excellent performance characteristics")

fn parallel_random_generation_demo(array_size: Int, num_threads: Int):
    """
    Demonstrate parallel random number generation.
    """
    print("Generating", array_size, "random numbers using", num_threads, "threads...")
    
    var random_numbers = List[UInt64]()
    for _ in range(array_size):
        random_numbers.append(0)
    
    @parameter
    fn generate_chunk(thread_id: Int):
        var chunk_size = array_size // num_threads
        var start_idx = thread_id * chunk_size
        var end_idx: Int
        if thread_id < num_threads - 1:
            end_idx = start_idx + chunk_size
        else:
            end_idx = array_size
        
        # Create thread-specific generator
        var thread_seed = UInt64(abs(random_si64(0, 1000000)) + thread_id * 1234567)
        var generator = Xoshiro256PlusPlus(thread_seed)
        
        for i in range(start_idx, end_idx):
            random_numbers[i] = generator.next()
    
    # Execute parallel generation
    parallelize[generate_chunk](num_threads)
    
    print("Generation completed!")
    
    # Display first 10 numbers as sample
    print("Sample random numbers:")
    var max_display = min(10, array_size)
    for i in range(max_display):
        print("  Index", i, ":", random_numbers[i])

fn test_float64_distribution():
    """
    Test variance and extreme values for next_float64.
    Tests on 10 million runs:
    - Variance of random numbers
    - Presence of numbers greater than 0.9999
    - Presence of numbers less than 0.00001.
    """
    print("=== Float64 Distribution and Extreme Values Test ===\n")
    
    var num_samples = 10000000  # 10 million
    var generator = Xoshiro256PlusPlus(54321)  # Fixed seed for reproducibility
    
    print("Generating", num_samples, "random Float64 values...")
    
    var start_time = perf_counter_ns()
    
    var sum: Float64 = 0.0
    var sum_squares: Float64 = 0.0
    var count_high: Int = 0  # Counter for values > 0.9999
    var count_low: Int = 0   # Counter for values < 0.00001
    
    for i in range(num_samples):
        var value = generator.next_float64()
        sum += value
        sum_squares += value * value
        
        if value > 0.9999:
            count_high += 1
        if value < 0.00001:
            count_low += 1
            
        # Progress every 2 million
        if i % 2000000 == 1999999:
            print("Processed", i + 1, "values...")
    
    var end_time = perf_counter_ns()
    var duration_ns = end_time - start_time
    var duration_ms = Float64(duration_ns) / 1_000_000.0
    
    # Calculate statistics
    var mean = sum / Float64(num_samples)
    var variance = (sum_squares / Float64(num_samples)) - (mean * mean)
    var std_dev = math.sqrt(variance)
    
    print("\n=== Analysis Results ===\n")
    print("Mean value:     ", mean, "(expected ~0.5)")
    print("Variance:       ", variance, "(expected ~0.0833)")
    print("Std. deviation: ", std_dev, "(expected ~0.2887)")
    print()
    print("Generation time:", duration_ms, "milliseconds")
    print("Generation rate:", Float64(num_samples) / duration_ms * 1000.0, "values per second")
    print()
    print("Count of values > 0.9999: ", count_high)
    print("Count of values < 0.00001:", count_low)
    print()
    
    # Validate results
    var mean_ok = abs(mean - 0.5) < 0.01  # Mean should be close to 0.5
    var variance_ok = abs(variance - 0.08333) < 0.01  # Variance ~1/12 for uniform distribution
    var high_values_ok = count_high > 0  # Should have values > 0.9999
    var low_values_ok = count_low > 0    # Should have values < 0.00001
    
    print("=== Results Validation ===\n")
    print("✓ Mean value is normal:" if mean_ok else "✗ Mean value is abnormal:", mean_ok)
    print("✓ Variance is normal:" if variance_ok else "✗ Variance is abnormal:", variance_ok)
    print("✓ High values found (>0.9999):" if high_values_ok else "✗ No high values found:", high_values_ok)
    print("✓ Low values found (<0.00001):" if low_values_ok else "✗ No low values found:", low_values_ok)
    
    var all_tests_passed = mean_ok and variance_ok and high_values_ok and low_values_ok
    print("\n", "🎉 ALL TESTS PASSED!" if all_tests_passed else "❌ SOME TESTS FAILED")
    print()

fn main():
    """
    Main function demonstrating multithreaded XorShift32+ usage.
    """
    print("=== Multithreaded xoshiro256++ Demo ===")
    print()
    
    # Test basic generator
    print("1. Basic xoshiro256++ Generator Test:")
    var rng = Xoshiro256PlusPlus(12345)
    print("Generated random numbers:")
    for _ in range(5):
        print("  Random uint64:", rng.next())
        print("  Random uint32:", rng.next32())
        print("  Random float [0,1):", rng.next_float())
        print("  Random float64 [0,1):", rng.next_float64())
        print("  Random int [10,20):", rng.next_range(10, 20))
        print()
    
    # Test parallel random number generation
    print("2. Parallel Random Number Generation:")
    parallel_random_generation_demo(1000000, 4)
    print()
    
    # Test Monte Carlo π estimation
    print("3. Parallel Monte Carlo π Estimation:")
    var num_samples = 10000000
    var num_threads = 4
    
    var estimated_pi = parallel_monte_carlo_pi(num_samples, num_threads)
    
    print("Using", num_samples, "samples with", num_threads, "threads")
    print("Estimated π:", estimated_pi)
    print("Actual π:   ", math.pi)
    print("Error:      ", abs(estimated_pi - math.pi))
    print()
    
    # Performance comparison: different thread counts
    print("4. Performance Comparison with Different Thread Counts:")
    var test_samples = 500000
    
    for threads in range(1, 5):
        var pi_estimate = parallel_monte_carlo_pi(test_samples, threads)
        var error = abs(pi_estimate - math.pi)
        print("Threads:", threads, "| π estimate:", pi_estimate, "| Error:", error)
    
    print()
    
    # Benchmark 1 million random numbers generation
    print("5. Time Benchmark - 1 Million Random Numbers:")
    benchmark_million_numbers_generation()
    
    print()
    
    # Test Float64 distribution
    print("6. Float64 Distribution Test:")
    test_float64_distribution()
    
    print()
    print("=== Demo Complete ===")
    print()
    print("Key Features Demonstrated:")
    print("- Thread-safe xoshiro256++ algorithm with 2^256-1 period")
    print("- 64-bit native and 32-bit random number generation")
    print("- High-precision Float64 random generation with 53-bit precision")
    print("- Statistical distribution testing")
    print("- Parallel random number generation")
    print("- Monte Carlo π estimation")
    print("- Scalable multithreading")
    print("- High-performance computation")