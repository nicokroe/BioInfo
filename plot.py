from pathlib import Path
import csv
import matplotlib.pyplot as plt


def create_plot(gProt: dict[str, tuple[float, float]], kegg_pathways: dict[str, list[str]]):
    """Creates a scatter plot. The function searches for kegg findings in gProt.

    Args:
        gProt (dict[str, tuple[float, float]]): A dict containing all proteins and their values
        kegg_pathways (dict[str, list[str]]): A dict containing annotated clusters of proteins
    """
    plots: dict[str, tuple[list[float], list[float]]] = {}
    for kegg, data in kegg_pathways.items():
        xs: list[float] = []
        ys: list[float] = []
        for protein in data:
            x, y = gProt[protein]
            xs.append(x)
            ys.append(y)
        plots[kegg] = (xs, ys)

    plt.xlabel("log2 Ratio H/L")
    plt.ylabel("log10 (protein abundance)")

    for kegg, data in plots.items():
        plt.scatter(data[0], data[1], label=kegg)

    plt.legend(loc="upper right")
    plt.show()


def create_kegg_pathways(read_tsv):
    # Read in rows from kegg data
    kegg_pathways: dict[str, list[str]] = {}
    for row in read_tsv:
        proteins = row[-1].split(",")
        kegg_pathways[row[0]] = proteins
    # Remove table headers from dict
    del kegg_pathways["#term ID"]
    return kegg_pathways


def create_gprot_dict(dataset2):
    # Read in rows from gProt
    gProt: dict[str, tuple[float, float]] = {}
    with dataset2.open() as file:
        # Skip headers of table
        next(file)
        for line in file:
            line = line.strip().split()
            name, input_value, property_value = line[1], line[3], line[4]
            gProt[name] = (float(input_value), float(property_value))
    return gProt


def main():
    dataset1 = Path(__file__).parent / "kegg.tsv"
    dataset2 = Path(__file__).parent / "log10.txt"

    tsv_reader = csv.reader(dataset1.open(), delimiter="\t")
    kegg_pathways = create_kegg_pathways(tsv_reader)
    gProt = create_gprot_dict(dataset2)

    create_plot(gProt, kegg_pathways)


if __name__ == "__main__":
    main()
