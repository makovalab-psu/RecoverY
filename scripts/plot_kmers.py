import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os

def plot_kmers() :

    plt.rcParams['figure.figsize'] = 12, 8
    plt.rcParams['xtick.labelsize'] = 14
    plt.rcParams['ytick.labelsize'] = 14

    ip_dir = "data"
    kmers_from_reads = ip_dir + "/kmers_from_reads"
    trusted_DSK_counts_acc_to_fsY = ip_dir + "/trusted_DSK_counts_acc_to_fsY"

    op_dir = "output"
    op_file = op_dir + "/kmerplot.png"

    #print "Generating k-mer plot"

    # if folder "output" doesn't exist, create the folder
    if not os.path.exists(op_dir):
        os.makedirs(op_dir)

    # if file op_r1.fastq exists, remove it
    if os.path.isfile(op_file):
        os.remove(op_file)

    # plotting kmer abundances

    abundance_list_fsY = []
    with open(kmers_from_reads, "r") as f:
        for line in f:
            abundance = float(str(line).strip().split()[1])
            abundance_list_fsY.append(abundance)

    abundance_list_trusted_genes = []
    with open(trusted_DSK_counts_acc_to_fsY, "r") as f:
        for line in f:
            abundance = float(str(line).strip().split()[1])
            abundance_list_trusted_genes.append(abundance)

    sns.set_style("white")

    bins_master = np.linspace(0, 1000, 1000)

    (n1, bins1, patches1) = plt.hist(abundance_list_trusted_genes, bins=bins_master, label='hst')

    (n2, bins2, patches2) = plt.hist(abundance_list_fsY, bins=bins_master, label='hst')

    plt.clf()

    hum_gene_kmers_counts = []
    with open(trusted_DSK_counts_acc_to_fsY, "r") as f:
        for line in f:
            abundance = float(str(line).split(' ')[1])
            if abundance > 0.0:
                hum_gene_kmers_counts.append(abundance)

    h = np.array(hum_gene_kmers_counts)
    five_pc = np.percentile(h, 5)
    ninety_five_pc = np.percentile(h, 95)

    xs = [i for i in range(0, 999, 1)]
    len(xs)

    sns.set_style("white")

    plt.plot(xs, n2, sns.xkcd_rgb["denim blue"], lw=5, label="K-mers from enriched data")
    plt.plot(xs, n1, sns.xkcd_rgb["pale red"], lw=5, label="Trusted-gene-kmers", )

    plt.legend(loc='upper right', fontsize=24)
    axes = plt.gca()

    # increase size of labels
    plt.tick_params(labelsize=24)

    # add commas
    axes.get_yaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))
    axes.get_xaxis().set_major_formatter(plt.FuncFormatter(lambda x, loc: "{:,}".format(int(x))))

    # set limits
    y_upper = 100000
    axes.set_xlim([3, ninety_five_pc + 50])
    axes.set_ylim([0, y_upper])

    # draw a vline
    plt.plot((five_pc, five_pc), (0, y_upper), sns.xkcd_rgb["black"], lw=4, linestyle='--', label="Trusted-genes")

    plt.xlabel('Abundance', fontsize=24)
    plt.ylabel('Number of kmers', fontsize=24)

    # adding a panel label
    # axes.text(-0.16, 1.03, "B", transform=axes.transAxes,fontsize=32, fontweight='bold', va='top', ha='right')


    # plt.show()
    plt.savefig(op_file, bbox_inches='tight')
