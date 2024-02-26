import os
import json
import tabulate


def pretty_time_diff_from_ns(nanoseconds):
    units = [
        (1, "ns"),  # nanoseconds
        (1e3, "μs"),  # microseconds
        (1e6, "ms"),  # milliseconds
        (1e9, "s"),  # seconds
        (60e9, "min"),  # minutes
        (3600e9, "h"),  # hours
        (86400e9, "d")  # days
    ]

    for factor, unit in units:
        if nanoseconds < factor:
            break
    time_amount = nanoseconds / (factor / 1000)
    return f"{time_amount:.2f} {unit}"


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
    # print(comparisons)
    headers = ["command", "bedtools time", "granges time", "granges speedup (%)"]
    table = []
    for comparison, tools in comparisons.items():
        if len(tools) == 2:
            bedtools_time, granges_time = tools.values()
            # Calculate the ratio and convert it to a percentage
            percent_faster = ((tools["bedtools"] - tools["granges"]) / tools["bedtools"]) * 100
            # Format the output to show 3 decimal places
            # print(f"{comparison} - GRanges is {percent_faster:.3f}% faster than Bedtools")
            table.append([comparison, pretty_time_diff_from_ns(tools["bedtools"]), 
                          pretty_time_diff_from_ns(tools["granges"]), percent_faster])
    print(tabulate.tabulate(table, headers=headers))


def main():
    criterion_dir = "target/criterion"
    calculate_ratios(criterion_dir)


if __name__ == "__main__":
    main()
