import os
import sys
import subprocess
import time
import random

def generate_points(num_points, filename):
    """Generates a file with a specified number of random 3D points."""
    with open(filename, 'w') as f:
        for _ in range(num_points):
            x = random.uniform(0, 1)
            y = random.uniform(0, 1)
            z = random.uniform(0, 1)
            f.write(f"{x} {y} {z}\n")

def main():
    voro_executable = sys.argv[1] if len(sys.argv) >= 2 else  './build/Voro'
    if not os.path.exists(voro_executable):
        print(f"Error: Executable not found at '{voro_executable}'")
        sys.exit(1)

    point_counts = [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]
    
    print("--- Starting Voro Benchmark ---")
    
    for n in point_counts:
        print(f"\nBenchmarking with {n:,} points...")
        
        points_filename = f"points_{n}.txt"
        
        generate_points(n, points_filename)

        # 1. Run voro++
        print(f"  Running voro++...")
        start_time_voropp = time.time()
        try:
            subprocess.run(
                ['voro++', '-g', '-c', '%i %x %y %z %n %V', '0', '1', '0', '1', '0', '1', points_filename],
                capture_output=True, text=True, check=True
            )
            elapsed_time_voropp = time.time() - start_time_voropp
            print(f"  voro++ took: {elapsed_time_voropp:.4f} seconds")
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"  Error running voro++: {e}")
            os.remove(points_filename)
            continue

        # 2. Run Voro tessellation
        print(f"  Running Voro tessellation...")
        start_time_voro = time.time()
        try:
            process = subprocess.run(
                [voro_executable, points_filename],
                capture_output=True,
                text=True,
                check=True
            )
            end_time_voro = time.time()
            
            elapsed_time_voro = end_time_voro - start_time_voro
            print(f"  Voro tessellation took: {elapsed_time_voro:.4f} seconds")
            # print("  Voro output:\n" + process.stdout)
            
        except subprocess.CalledProcessError as e:
            print(f"  Error running Voro for {n} points.")
            print(f"  Return code: {e.returncode}")
            print(f"  Stdout: {e.stdout}")
            print(f"  Stderr: {e.stderr}")
            # Stop benchmark if one run fails
            break
        except FileNotFoundError:
            print(f"  Error: The command '{voro_executable}' was not found.")
            break

        # 3. Clean up the points file
        # print(f"  Cleaning up '{points_filename}'...")
        os.remove(points_filename)

    print("\n--- Benchmark Finished ---")

if __name__ == "__main__":
    main()
