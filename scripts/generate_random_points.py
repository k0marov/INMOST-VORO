#!/usr/bin/env python3
import argparse
import random
import sys

def generate_random_points(num_points, filename):
    print(f"Generating {num_points} random points and saving to '{filename}'...")
    written = 0

    try:
        with open(filename, 'w') as f:
            for i in range(num_points):
                x = random.uniform(0, 1)
                y = random.uniform(0, 1)
                z = random.uniform(0, 1)
                f.write(f"{x} {y} {z}\n")
                written += 1

                if written % 100000 == 0:
                    print(f"  {written}/{num_points} points written.")
    except IOError as e:
        print(f"Error writing to file: {e}", file=sys.stderr)
        sys.exit(1)

    print("Generation complete.")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate a file containing random 3D points in [0,1]."
    )
    parser.add_argument(
        "-n", "--num-points",
        type=int,
        required=True,
        help="Number of points to generate."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output filename."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.num_points <= 0:
        print("Error: --num-points must be a positive integer.", file=sys.stderr)
        sys.exit(1)

    generate_random_points(args.num_points, args.output)
