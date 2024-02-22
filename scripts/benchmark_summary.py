import os
import json


def extract_mean_point_estimate(estimates_path):
    with open(estimates_path, "r") as file:
        data = json.load(file)
    return data["mean"]["point_estimate"]


def calculate_ratios(criterion_dir):
    comparisons = {}  # Store comparison: ratio

    subcommand_benches = [d for d in os.listdir(criterion_dir) if d != "report"]

    for subcommand in subcommand_benches:
        benches = dict()
        runs = os.listdir(os.path.join(criterion_dir, subcommand))

        for dir in runs:
            if dir.startswith("report"):
                continue
            if dir.startswith("bedtools"):
                bedtools_bench = os.path.join(
                    criterion_dir, subcommand, dir, "new", "estimates.json"
                )
                benches["bedtools"] = extract_mean_point_estimate(bedtools_bench)
            if dir.startswith("granges"):
                granges_bench = os.path.join(
                    criterion_dir, subcommand, dir, "new", "estimates.json"
                )
                benches["granges"] = extract_mean_point_estimate(granges_bench)
            comparisons[subcommand] = benches

    # Calculate and print ratios for each comparison
    for comparison, tools in comparisons.items():
        if len(tools) == 2:
            bedtools_time, granges_time = tools.values()
            # Calculate the ratio and convert it to a percentage
            percent_faster = ((bedtools_time - granges_time) / bedtools_time) * 100
            # Format the output to show 3 decimal places
            print(
                f"{comparison} - GRanges is {percent_faster:.3f}% faster than Bedtools"
            )


def main():
    criterion_dir = "target/criterion"
    calculate_ratios(criterion_dir)


if __name__ == "__main__":
    main()
