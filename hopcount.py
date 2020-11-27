#! /usr/bin/env python3

"""
Python script for the analysis of transposition events in sequencing data

INPUT:
    - A reference GenBank file, with genic annotations. Designed with genbank files from NCBI
    - A folder, at ./fastq/, that contains the trimmed reads ready for processing in .fastq or .fastq.gz format
    - The configuration file, at ./hopcount.conf, that contains the processing options for hopcount.py
OUTPUT:
    - Each output file will have a prefix derived from the read file, in its own folder in ./output/
    - _hopcount.tsv, which is a record of each transposition event recorded along with its determined genetic context
    - _aggregate.tsv, which aggregates the records in hopcount.tsv into each region entry, with DvalGenome calculations and more
    - _sense.tsv. which analyzes contributions in DvalGenome regarding same-sense insertions between adjacent genes
    - _antisense.tsv, which analyzes contributions in DvalGenome regarding opposite-sense insertions between adjacent genes
    - .FASTA and .GFF3 files from the reference
    - Optionally, the BAM or SAM file from bowtie2 mapping
    - Optionally, WIG files for visualization into a genome browser

REQUIREMENTS:
    - Python 3.6+
    - BioPython
    - bcbio-gff 0.6.6 (https://pypi.org/project/bcbio-gff/)

Author: jacob.bourgeois@tufts.edu
Affiliation - Camilli Laboratory, Tufts University Graduate School of Biomedical Sciences, Molecular Microbiology

Version: 1.0
- Rewrite of the flow of the program to be able automatically bin reads that fail percent thresholding to adjacent
 intergenic regions
- Fixes various bugs in reporting
- Various optimizations
- Now reports approximate complexity, defined as the number of CDS inserted over the total possible CDS, and other stats
- This is intended to be the first stable release to accompany the CSHL Handbook
"""


import csv
from collections import defaultdict
from statistics import mean
import subprocess
import os
import math

# Check if Biopython is installed
try:
    from Bio import SeqIO
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("Error! Biopython is not installed. Please install Biopython [https://biopython.org/wiki/Download]")
    quit()

# Check if experimental GFF module is installed
try:
    from BCBio import GFF
except ImportError:
    print("Error! GFF module is not installed. Please install bcbio gff tools [https://pypi.org/project/bcbio-gff/]")
    quit()


# use samtools to convert sam to bam
def bamify(sam, bam):
    subprocess.run(['samtools', 'view', '-b', sam, '-o', bam, '-@', '4'], check=True)
    os.remove(sam)  # delete old sam file
    return


# converts gbk to fasta for bowtie
def gbk_to_fasta(gbk, fasta):

    with open(gbk, 'r') as handle:
        records = list(SeqIO.parse(handle, 'gb'))
        with open(fasta, 'w') as out_handle:
            SeqIO.write(records, out_handle, "fasta")
    return


# uses bowtie2-build to make a reference index
def bowtie2_build(ref, ind):
    subprocess.run(['bowtie2-build', ref, ind], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    return


# align using bowtie2
# unpaired mode
# note bowtie2 uses stderr for output, oddly enough
def bowtie2_unpaired(ind, fq, sam, ops):

    # Very rarely, there is a poorly formatted FASTQ file that catches. Return a fail.

    # Convert options dictionary to a list
    command = ['bowtie2', '-x', ind, '-U', fq, '-S', sam]
    command.extend([item for k in ops for item in (k, ops[k])])

    try:
        proc = subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
                              universal_newlines=True)
        print('\n', proc.stderr)  # Display the mapping statistics in stderr.
        return proc.stderr
    except subprocess.CalledProcessError as e:
        print("Bowtie2 error! {0}".format(e))
        print(e.stderr)
        return False


# Get CDS info from a genbank file
def extract_gene_info(gbk, is_circular='true'):

    # Open the gbk file
    gbk_prefix = os.path.splitext(os.path.basename(gbk))[0]

    with open(gbk, 'r') as input_handle:

        # We do not know ahead of time how many elements (ie. chromosomes) are present.
        # So, define a dict to return, indexed by chromosome name, like the gbk file
        genbank_elements = dict()
        genome_lengths = dict()
        gff_records = []
        split_locus_tags = dict()
        cds_length_data = dict()

        # Foe each chromosome in the input GBK file
        for record in SeqIO.parse(input_handle, 'gb'):

            # grab genbank info
            chr_name = record.id
            description = record.description
            cds = [k for k in record.features if k.type == 'CDS']  # Standard convention from NCBI GBK
            cds_length_data[chr_name] = get_region_lengths(cds)

            total_length = len(record.seq)

            print("Extracted {0} CDS from {1} ({2})".format(len(cds), chr_name, description))

            # If a genome segment has zero coding sequences and we are circular,
            # then print a warning and add start/end keys (v0.53)
            if len(cds) == 0 and is_circular == 'true':
                print("Warning! Element {0} has no coding sequences on circular chromosome. Creating start/end keys...".format(chr_name), end='')
                start_of_chromosome = SeqFeature(FeatureLocation(0, 0), qualifiers={'locus_tag': ['Start of chromosome {0}'.format(chr_name)],
                                                                          'gene': ['Start of chromosome {0}'.format(chr_name)], 'product': ['N/A'],
                                                                          'protein_id': ['N/A'], 'note': ['N/A']})
                end_of_chromosome = SeqFeature(FeatureLocation(total_length, total_length), qualifiers={'locus_tag': ['End of chromosome {0}'.format(chr_name)],
                                                                          'gene': ['End of chromosome {0}'.format(chr_name)], 'product': ['N/A'],
                                                                          'protein_id': ['N/A'], 'note': ['N/A']})
                cds = [start_of_chromosome, end_of_chromosome]
                print('done.')

            # Check for genes that break over the origin - rarely, some genbank files have a genic entry that splits
            # over the origin and screws up all the calculations
            multiple_part_genes = [k for k in cds if len(k.location.parts) > 1]
            for gene in multiple_part_genes:

                part1, part2 = gene.location.parts[0], gene.location.parts[1]
                if part1.start == 0 or part2.start == 0:

                    # Orient so part1 extends from the start of the chromosome
                    if part1.start != 0:
                        foo = part2
                        part2 = part1
                        part1 = foo

                        # The first part of the gene "starts" at zero if over the origin
                        print("Warning! Split origin gene detected for entry {0}".format(gene.qualifiers['locus_tag'][0]))

                        # Extract end coordinate for the current gene and the start coordinate of the next gene
                        # We have to make this a compound location entry

                        part1_gene_start = 0
                        part1_gene_end = part1.end
                        part2_gene_start = part2.start
                        part2_gene_end = total_length

                        feature_location_part1 = FeatureLocation(part1_gene_start, part1_gene_end)
                        feature_location_part2 = FeatureLocation(part2_gene_start, part2_gene_end)

                        start_entry = SeqFeature(feature_location_part1, type='CDS',
                                                qualifiers={k: gene.qualifiers[k] for k in gene.qualifiers})
                        end_entry = SeqFeature(feature_location_part2, type='CDS',
                                                qualifiers={k: gene.qualifiers[k] for k in gene.qualifiers})
                        cds.remove(gene)
                        cds.insert(0, start_entry)
                        cds.append(end_entry)

            # Assign CDS and chromosome information to the master CDS dictionary
            genbank_elements[chr_name] = cds
            genome_lengths[chr_name] = {'TotalLength': total_length}

            # If the genome isn't circular, then we need to anchor reads that fall before the start CDS or after the
            # end CDS with an artificial key.
            if is_circular == 'false':
                print('Creating start and end chromosomal keys...')
                start_of_chromosome = SeqFeature(FeatureLocation(0, 0), qualifiers={'locus_tag': ['Start of chromosome {0}'.format(chr_name)],
                                                                          'gene': ['Start of chromosome {0}'.format(chr_name)], 'product': ['N/A'],
                                                                          'protein_id': ['N/A'], 'note': ['N/A']})
                end_of_chromosome = SeqFeature(FeatureLocation(total_length, total_length), qualifiers={'locus_tag': ['End of chromosome {0}'.format(chr_name)],
                                                                          'gene': ['End of chromosome {0}'.format(chr_name)], 'product': ['N/A'],
                                                                          'protein_id': ['N/A'], 'note': ['N/A']})

                genbank_elements[chr_name].insert(0, start_of_chromosome)
                genbank_elements[chr_name].append(end_of_chromosome)

            # Add CDS element information to GFF records, needed to make the GFF3 file for browser viewing
            gff_records.append(SeqRecord(record.seq, chr_name, chr_name, description, features=cds))

    # Write a GFF3 file for data visualization in a genome browser
    # Uses an experimental module from BioPython, so this command may become deprecated
    with open('{0}.gff3'.format(gbk_prefix), 'w') as handle:
        GFF.write(gff_records, handle)
        print("Wrote GFF3 file to {0}".format(os.path.join(os.getcwd(), '{0}.gff3'.format(gbk_prefix))))

    return genbank_elements, genome_lengths, split_locus_tags, cds_length_data


# Create a start, end, strand, length dictionary
# This is mainly for handling CDS regions split over the origin
def get_region_lengths(genes):

    size_dict = dict()

    for gene in genes:

        # get start, end, strand, length
        start = gene.location.start
        end = gene.location.end
        length = end - start
        strand = gene.location.strand
        locus = gene.qualifiers['locus_tag'][0]

        # Do a fix for split origin CDS
        if len(gene.location.parts) > 1:
            if gene.location.parts[0].start == 0:
                print("Fixing split CDS {0}...".format(locus))
                start = gene.location.parts[1].start
                end = gene.location.parts[0].end
                length = gene.location.parts[0].end + (gene.location.parts[1].end - gene.location.parts[1].start)

        size_dict[locus] = {'Start': start, 'End': end, 'Length': length, 'Strand': strand}
    return size_dict


# Function that generates a dictionary of all possible IG regions based on CDS data
def generate_ig_regions(cds_data, is_circular, genome_lengths):

    ig_region_lengths = dict()
    for chromosome in cds_data:
        # print("Generating keys for {0}...".format(chromosome))

        # holder for inserted entries
        ig_entries = dict()

        # iterate through the list of genes in a given chromosome
        total_genes = len(cds_data[chromosome])
        for k in range(0, total_genes - 1):
            current_gene = cds_data[chromosome][k]
            next_gene = cds_data[chromosome][k + 1]

            # Extract end coordinate for the current gene and the start coordinate of the next gene
            current_gene_end = current_gene.location.end
            next_gene_start = next_gene.location.start

            if next_gene_start > current_gene_end:
                # There is space, thus, an intergenic region!

                locus = '{0}--{1}_ig'.format(current_gene.qualifiers['locus_tag'][0],
                                             next_gene.qualifiers['locus_tag'][0])

                product = 'N/A'

                protein_id = 'N/A'

                gene = ''

                note = 'Intergenic region between {0} and {1}'.format(current_gene.qualifiers['locus_tag'][0],
                                                                      next_gene.qualifiers['locus_tag'][0])

                ig_entry = SeqFeature(FeatureLocation(current_gene_end, next_gene_start, strand=None), type='IG',
                                      qualifiers={'locus_tag': [locus],
                                                  'gene': [gene],
                                                  'product': [product],
                                                  'protein_id': [protein_id],
                                                  'note': [note]})

                ig_entries[ig_entry] = k + 1

        # Now, let's insert the ig entries into the master list
        # Remember, everytime we add an entry, the index shifts up by one!
        shift = 0
        for entry in ig_entries:
            cds_data[chromosome].insert(ig_entries[entry] + shift, entry)
            shift += 1
        # print("Added {0} intergenic entries to {1}".format(len(ig_entries), chromosome))
        ig_region_lengths[chromosome] = get_region_lengths(ig_entries)

        # Technically there is a final possible intergenic region between the end of the final CDS and start of the
        # next CDS if the genome is circular, and we aren't at the end of the chromosome (by an existing split locus)
        # We can detect this case if the last CDS end is not the length of the chromosome

        if is_circular == 'true' and cds_data[chromosome][len(cds_data[chromosome]) - 1].location.end < genome_lengths[chromosome]['TotalLength']:
            # print("There is room at the end for a split entry!")

            # Extract end coordinate for the current gene and the start coordinate of the next gene
            # We have to make this a compound location entry

            start_gene = cds_data[chromosome][0]
            end_gene = cds_data[chromosome][len(cds_data[chromosome]) - 1]

            part1_gene_start = 0
            part1_gene_end = start_gene.location.start
            part2_gene_start = end_gene.location.end
            part2_gene_end = genome_lengths[chromosome]['TotalLength']

            locus = '{0}--{1}_ig'.format(end_gene.qualifiers['locus_tag'][0],
                                         start_gene.qualifiers['locus_tag'][0])
            product = 'N/A'
            protein_id = 'N/A'
            gene = ''
            note = 'Intergenic region between {0} and {1}'.format(end_gene.qualifiers['locus_tag'][0],
                                                                  start_gene.qualifiers['locus_tag'][0])

            feature_location_part1 = FeatureLocation(part1_gene_start, part1_gene_end)
            feature_location_part2 = FeatureLocation(part2_gene_start, part2_gene_end)

            ig_entry_1 = SeqFeature(feature_location_part1, type='IG',
                                  qualifiers={'locus_tag': [locus],
                                              'gene': [gene],
                                              'product': [product],
                                              'protein_id': [protein_id],
                                              'note': [note]})
            ig_entry_2 = SeqFeature(feature_location_part2, type='IG',
                                    qualifiers={'locus_tag': [locus],
                                                'gene': [gene],
                                                'product': [product],
                                                'protein_id': [protein_id],
                                                'note': [note]})
            cds_data[chromosome].insert(0, ig_entry_1)
            cds_data[chromosome].append(ig_entry_2)

            size_data_payload = {'Start': part2_gene_start, 'End': part1_gene_end, 'Strand': None, 'Length': part1_gene_end + (part2_gene_end - part2_gene_start)}
            ig_region_lengths[chromosome][locus] = size_data_payload

    return cds_data, ig_region_lengths


# Obtains a transposition read histogram from a SAM file
def get_transp_read_pos(sam_file, qual_cutoff):

    # Nasty layered defaultdict to handle multiple genomic elements in one SAM file
    pos_histogram = defaultdict(lambda: defaultdict(lambda: [0, 0, []]))  # Position histogram of transposition events.

    rows_included = 0
    total_rows = 0
    match_rows_included = 0

    with open(sam_file, 'r') as reads:

        # Or, maybe there's a way to iterate over stdout buffer with samtools.

        for read in reads:

            # shed the header lines
            while read[0] == '@':
                read = next(reads)

            # First, split the read line by tab
            read = read.split('\t')
            total_rows += 1

            # Now, extract the read position and direction. We don't have column headers, but by default...
            read_strand_bitval = read[1]
            read_chromosome = read[2]
            read_pos = read[3]
            read_seq = read[9]
            read_length = len(read_seq)
            read_mapq = int(read[4])

            # Do quality cutoff

            if read_mapq >= qual_cutoff:
                rows_included += 1

                # Regarding the bitval, 0 is pos strand, and 16 is reverse strand, if the reads aligned
                if read_strand_bitval == '0':

                    # If the read_strand is positive, the read_pos obtained is proximal to the transposition site
                    # Thus, add one to the position in the read
                    # pos_histogram[read_pos][1] += 1
                    pos_histogram[read_chromosome][read_pos][0] += 1
                    pos_histogram[read_chromosome][read_pos][2].append(read_mapq)
                    match_rows_included += 1

                elif read_strand_bitval == '16':

                    # If the read strand is negative, then the read_pos is distal to the transposition site
                    # Thus, add the read length to the obtained read_pos, and index that in the dict
                    transp_read_pos = str(int(read_pos) + read_length)
                    pos_histogram[read_chromosome][transp_read_pos][1] += 1
                    pos_histogram[read_chromosome][transp_read_pos][2].append(read_mapq)
                    match_rows_included += 1

                else:
                    read_strand = 'UNKNOWN'
                    # print("Warning! Unexpected read flag: {0}".format(read_strand_bitval))

    # Now for each pos, go ahead and merge the fastq scores
    for chrom in pos_histogram:
        for entry in pos_histogram[chrom]:
            pos_histogram[chrom][entry][2] = mean(pos_histogram[chrom][entry][2])

    hop_stats = "{0} of {1} reads passed MAPQ quality threshold of {2} ({3:.2f}%)".format(rows_included, total_rows, qual_cutoff, 100*rows_included/total_rows)
    print(hop_stats)

    return pos_histogram, hop_stats


# Function determines genomic context of given position
#
# Lots of try-except clauses to catch NCBI errors
def get_genomic_context(genes, pos):

    # We have to set up context as a defaultdict because there may exist transposition events in which gene annotations
    # overlap. For these, we want to add one to both entries.
    context = defaultdict(dict)

    # first, let's get the CDS/IG regions in which these overlap using list comprehension
    gene_hits = [k for k in genes if (k.location.start <= pos <= k.location.end)]

    if gene_hits:
        i = 0
        for cds in gene_hits:

            i += 1
            try:
                locus = cds.qualifiers['locus_tag'][0]

            except KeyError:
                locus = 'N/A'

            try:
                genei = cds.qualifiers['gene'][0]
            except KeyError:
                genei = ''

            try:
                product = cds.qualifiers['product'][0]
            except KeyError:
                product = 'N/A'

            try:
                protein_id = cds.qualifiers['protein_id'][0]
            except KeyError:
                protein_id = 'N/A'

            try:
                note = cds.qualifiers['note'][0]
            except KeyError:
                note = 'N/A'

            context['gene{0}'.format(i)] = {'locus': locus,
                                            'gene': genei,
                                            'product': product,
                                            'protein_id': protein_id,
                                            'note': note,
                                            'END_LOCUS_POS': cds.location.end}

    return context


# Returns a new gene contextual list that includes percent thresholding parameters
def threshold_gene_list(input_gene_list, percent_threshold, is_circular, genome_lengths):

    thresholded_entries = dict()

    for chromosome in input_gene_list:

        thresholded_entries[chromosome] = []

        for entry in input_gene_list[chromosome]:

            start = entry.location.start
            end = entry.location.end
            threshold_start = start + math.ceil(percent_threshold * (end - start))
            threshold_end = end - math.ceil(percent_threshold * (end - start))

            # Create a new SeqFeature
            new_entry = SeqFeature(FeatureLocation(threshold_start, threshold_end), type='CDS',
                                   strand=entry.strand, qualifiers={k: v for k, v in entry.qualifiers.items()})
            thresholded_entries[chromosome].append(new_entry)

        # generate a new intergenic library
        thresholded_entries, _ = generate_ig_regions(thresholded_entries, is_circular, genome_lengths)

    return thresholded_entries


# Function weaves together gene data and position histogram to make 'hopcount' file
# Perform percent-within thresholding here
def write_hopcount_file(out, thresholded_gene_list, histogram, original_gene_list, percent_cutoff, genome_lengths):

    with open(out, 'w') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['Reference', 'Position', 'Locus', 'Gene', 'PlusCount', 'MinusCount', 'Total', 'MAPQ', 'Product',
                        'ProteinID', 'Note'])

        total_positions = 0
        tossed_positions = 0
        write_positions = 0
        multiple_entries = 0
        total_unique_entries = set()
        total_genome_length = 0

        for chrom in sorted(histogram):

            original_cds = [k.qualifiers['locus_tag'][0] for k in original_gene_list[chrom]]
            total_positions += len(histogram[chrom])
            total_genome_length += int(genome_lengths[chrom]['TotalLength'])

            for pos in sorted(histogram[chrom], key=int):

                # The transposition histogram is sorted; therefore we climb our CDS list in order as well.

                pos_strand_events = histogram[chrom][pos][0]
                minus_strand_events = histogram[chrom][pos][1]
                total_strand_events = pos_strand_events + minus_strand_events

                gene_context_info = get_genomic_context(thresholded_gene_list[chrom], int(pos))
                if len(gene_context_info) > 1:
                    multiple_entries += 1

                for context in gene_context_info:

                    # Occasionally, for genes that are very close together, thresholding by percent within can create
                    # new intergenic regions that didn't originally exist. Therefore, we check to make sure they
                    # existed to begin with, and only write those to file.

                    if gene_context_info[context]['locus'] in original_cds:

                        writer.writerow([chrom,
                                         pos,
                                         gene_context_info[context]['locus'],
                                         gene_context_info[context]['gene'],
                                         pos_strand_events,
                                         minus_strand_events,
                                         total_strand_events,
                                         histogram[chrom][pos][2],
                                         gene_context_info[context]['product'],
                                         gene_context_info[context]['protein_id'],
                                         gene_context_info[context]['note']])

                        write_positions += 1
                        total_unique_entries.add(pos)

                    else:
                        # This position landed in an intergenic region that arose from aritfical CDS shortening
                        tossed_positions += 1
                        pass

    insertional_density = total_genome_length / len(total_unique_entries)

    stat_payload = ["Total unique positions in histogram: {0}".format(total_positions),
                    "Wrote {0} entries to file, tossed {2} entries [PERCENT_WITHIN_READ_CUTOFF = {1}%]".format(write_positions, percent_cutoff*100, tossed_positions),
                    "{0} positions fell under two or more CDS regions".format(multiple_entries),
                    "Total unique positions written: {0}".format(len(total_unique_entries)),
                    "Approximate positional complexity (before read count thresholding): 1 insertion per {0:.2f} bp".format(
                        insertional_density)]

    for stat in stat_payload:
        print(stat)

    return stat_payload


# Function reads a hopcount file, aggregates them, and finds deterministic values
def aggregate_hopcount(hop_file, genome_chars, gene_info, read_count_threshold=0):

    total_positions = 0
    total_reads = 0
    kept_positions = 0
    kept_reads = 0

    with open(hop_file, 'r') as input_handle:
        reader = csv.DictReader(input_handle, delimiter='\t')

        agg_site_data = defaultdict(lambda: [0, 0, 0, 0, 0])
        agg_data = dict()
        agg_pos_data = dict()
        genome_count = defaultdict(int)

        for row in reader:

            total_positions += 1

            chromosome = row['Reference']
            pos = row['Position']
            locus = row['Locus']
            gene = row['Gene']
            plus_count = int(row['PlusCount'])
            minus_count = int(row['MinusCount'])
            product = row['Product']
            protein_id = row['ProteinID']
            note = row['Note']

            total_reads += (plus_count + minus_count)

            # Threshold based on total reads here - need at least X number of reads in plus and minus
            if (plus_count + minus_count) <= read_count_threshold:
                # Failed threshold, leave
                pass
            else:
                agg_site_data[locus][0] += 1  # sites counter
                agg_site_data[locus][1] += plus_count  # plus read counter
                agg_site_data[locus][2] += minus_count  # minus read counter

                if locus not in agg_pos_data:
                    agg_pos_data[locus] = dict()
                    agg_pos_data[locus]['DATA'] = dict()

                agg_pos_data[locus]['DATA'][pos] = plus_count + minus_count
                agg_pos_data[locus]['SUM'] = sum(agg_pos_data[locus]['DATA'].values())
                agg_pos_data[locus]['MAX'] = max(agg_pos_data[locus]['DATA'], key=lambda k: agg_pos_data[locus]['DATA'][k])
                agg_pos_data[locus]['CONTRIB'] = {agg_pos_data[locus]['MAX']: '{0:.2f}'.format(100*agg_pos_data[locus]['DATA'][agg_pos_data[locus]['MAX']] / agg_pos_data[locus]['SUM'])}

                if plus_count >= 1:
                    agg_site_data[locus][3] += 1  # unique plus sites counter
                if minus_count >= 1:
                    agg_site_data[locus][4] += 1  # unique minus sites counter

                agg_data[locus] = {'Reference': chromosome,
                                   'Locus': locus,
                                   'Sites': agg_site_data[locus][0],
                                   'PlusCount': agg_site_data[locus][1],
                                   'MinusCount': agg_site_data[locus][2],
                                   'PlusSites': agg_site_data[locus][3],
                                   'MinusSites': agg_site_data[locus][4],
                                   'TotalCount': agg_site_data[locus][1] + agg_site_data[locus][2],
                                   'LargestHopContrib': '{0}% from hop position {1}'.format(agg_pos_data[locus]['CONTRIB'][agg_pos_data[locus]['MAX']], agg_pos_data[locus]['MAX']),
                                   'ProteinID': protein_id,
                                   'Gene': gene,
                                   'Product': product,
                                   'Note': note}
                genome_count[chromosome] += (plus_count + minus_count)
                kept_positions += 1
                kept_reads += (plus_count + minus_count)

    for entry in genome_count:
        genome_chars[entry].update({'TotalCounts': genome_count[entry]})

    agg_stat = "Of {0} reads covering {1} position entries, {2} reads covering {3} positions were kept [READ_THRESHOLD = {4}].".format(total_reads, total_positions, kept_reads, kept_positions, read_count_threshold)
    print(agg_stat)

    # Create null entries for entries that did not have transposition events
    print("Making null transposition entries...")
    agg_data = make_null_transposition_events(agg_data, gene_info=gene_info)

    return agg_data, agg_stat


def make_null_transposition_events(agg_data, gene_info):

    for chromosome in gene_info:

        null_entries = [k.qualifiers['locus_tag'][0] for k in gene_info[chromosome] if k.qualifiers['locus_tag'][0] not in agg_data]

        for entry in null_entries:

            cds_info = [k for k in gene_info[chromosome] if k.qualifiers['locus_tag'][0] == entry][0]

            product = cds_info.qualifiers['product'][0]
            try:
                protein_id = cds_info.qualifiers['protein_id'][0]
            except KeyError:
                protein_id = 'N/A'
            try:
                gene = cds_info.qualifiers['gene'][0]
            except KeyError:
                gene = ''

            agg_data[entry] = {'Reference': chromosome,
                          'Locus': entry,
                          'PlusSites': 0,
                          'MinusSites': 0,
                          'Sites': 0,
                          'PlusCount': 0,
                          'MinusCount': 0,
                          'TotalCount': 0,
                          'ProteinID': protein_id,
                          'Gene': gene,
                          'Product': product}

    return agg_data


# Function appends gene length and other goodies into the aggregate dictionary
def align_gene_info(agg, chr_lengths):

    # Let's add the positional information for cds first
    # Specifically, start, end, strand, and length

    for entry in agg:
        # Append values into aggregate container
        agg[entry].update(chr_lengths[agg[entry]['Reference']][entry])

    # remove start and end chromosome entries if they exist
    for chrom in chr_lengths:
        try:
            agg.pop('Start of chromosome {0}'.format(chrom))
            agg.pop('End of chromosome {0}'.format(chrom))
        except KeyError:
            pass

    return


# Function calculates expectation values for aggregate data
# Also calculates a plus/minus ratio
def expect_val(agg, chr_lengths):

    # define subfunction for fitness calculations
    def expect(a, b, c, d): return (a/b) / (c/d)

    # use global values for chromosome attributes - what is the total length and counts across all chromosomes?
    absolute_length, absolute_counts = 0, 0
    for entry in chr_lengths:
        genome_total_count = chr_lengths[entry]['TotalCounts']
        genome_length = chr_lengths[entry]['TotalLength']
        absolute_counts += genome_total_count
        absolute_length += genome_length

    for entry in agg:

        # chromosome = agg[entry]['Reference']

        region_plus = agg[entry]['PlusCount']
        region_total_count = agg[entry]['TotalCount']

        # Region ratio
        if region_total_count > 0:
            agg[entry]['percentPlus'] = '{0:.3f}'.format(100 * region_plus / region_total_count)
        else:
            agg[entry]['percentPlus'] = ''

        # genome_total_count = chr_lengths[chromosome]['TotalCounts']
        region_length = agg[entry]['Length']
        # genome_length = chr_lengths[chromosome]['TotalLength']

        fitness = expect(region_total_count, absolute_counts, region_length, absolute_length)

        test = 1/(expect(1, absolute_counts, region_length, absolute_length))

        agg[entry]['DvalGenome'] = fitness
        agg[entry]['ExpectedCount'] = test

    return


# Function looks for antagonistic and synergistic relationships between intergenic insertions and downstream Dval
# Desperately in need of more efficient coding. Working at the moment.
def agg_analyze_combos(agg_data, gene_data, genome_length, is_circular):

    sense_data, antisense_data = dict(), dict()

    # Obtain a list of all applicable loci, sorted by position
    # Separate into CDS and IG regions

    for chrom in gene_data:

        # sorting in this way preserves the actual order of CDS by nt positions

        agg_chrom_data = {k: v for (k, v) in agg_data.items() if agg_data[k]['Reference'] == chrom}

        sorted_regions = sorted(agg_chrom_data.items(), key=lambda x: x[1]['Start'])
        loci = [k[0] for k in sorted_regions if '_ig' not in k[0]]
        ig = [k[0] for k in sorted_regions if k[0] not in loci]
        all_loci = loci + ig
        i = 0

        # for each entry
        for locus in all_loci:

            # Obtain the upstream and downstream gene data. In the case of intergenic regions, it is contained in the gene
            # name - use a string split operation
            # Use a try block in case of some weird end of chromosome shit

            upstream_gene, downstream_gene = '', ''

            if '_ig' in locus:
                upstream_gene, downstream_gene = locus.split('--')[0], locus.split('--')[1][:-3]
            else:
                try:
                    if i == 0:  # we are at the start of the chromosome.
                        # later, consider if genome is circular...the upstream gene is the last gene!
                        downstream_gene = loci[i+1]
                        if is_circular == 'true':
                            upstream_gene = loci[len(loci)-1]
                    elif i == len(loci) - 1:  # we are at the end of the chromosome
                        # if chromosome is circular, the downstream gene is the first gene!
                        upstream_gene = loci[i-1]
                        if is_circular == 'true':
                            downstream_gene = loci[0]
                    else:
                        upstream_gene = loci[i-1]
                        downstream_gene = loci[i+1]
                except KeyError:
                    print("KeyError! {0} does not have a vaild upstream or downstream genomic site!".format(locus))
                i += 1  # advance the list of CDS

            #print("{0}: up:{1}, down:{2}".format(locus, upstream_gene, downstream_gene))

            # get current gene data

            current_gene_reference = agg_data[locus]['Reference']
            current_gene_name = agg_data[locus]['Gene']
            current_gene_strand = agg_data[locus]['Strand']
            current_gene_plus_sites = agg_data[locus]['PlusSites']
            current_gene_minus_sites = agg_data[locus]['MinusSites']
            current_gene_total_sites = agg_data[locus]['Sites']
            current_gene_plus_count = agg_data[locus]['PlusCount']
            current_gene_minus_count = agg_data[locus]['MinusCount']
            current_gene_total_count = agg_data[locus]['TotalCount']
            current_gene_start = agg_data[locus]['Start']
            current_gene_end = agg_data[locus]['End']
            current_gene_length = agg_data[locus]['Length']
            current_gene_dval = agg_data[locus]['DvalGenome']
            current_gene_product = agg_data[locus]['Product']

            # Entries which had no hopes will have zero counts in total. Use ZeroDivisionError to escape.

            try:
                current_gene_percent_plus_count = 100*current_gene_plus_count/current_gene_total_count
                current_gene_plus_dval = current_gene_dval * (current_gene_percent_plus_count/100)
            except ZeroDivisionError:
                current_gene_percent_plus_count = 0
                current_gene_plus_dval = 0

            try:
                current_gene_percent_minus_count = 100*current_gene_minus_count/current_gene_total_count
                current_gene_minus_dval = current_gene_dval * (current_gene_percent_minus_count/100)
            except ZeroDivisionError:
                current_gene_percent_minus_count = 0
                current_gene_minus_dval = 0

            # Collect upstream and downstream gene info
            if upstream_gene and 'Start of chromosome' not in upstream_gene:
                upstream_gene_reference = agg_data[upstream_gene]['Reference']
                upstream_gene_name = agg_data[upstream_gene]['Gene']
                upstream_gene_strand = agg_data[upstream_gene]['Strand']
                upstream_gene_plus_sites = agg_data[upstream_gene]['PlusSites']
                upstream_gene_minus_sites = agg_data[upstream_gene]['MinusSites']
                upstream_gene_plus_count = agg_data[upstream_gene]['PlusCount']
                upstream_gene_minus_count = agg_data[upstream_gene]['MinusCount']
                upstream_gene_start = agg_data[upstream_gene]['Start']
                upstream_gene_end = agg_data[upstream_gene]['End']
                upstream_gene_length = agg_data[upstream_gene]['Length']
                upstream_gene_dval = agg_data[upstream_gene]['DvalGenome']
                upstream_gene_product = agg_data[upstream_gene]['Product']

                # intergenic distance hack
                if '_ig' in locus:
                    genic_distance = agg_data[locus]['Length']
                else:
                    # In case we are wrapping around the chromosome, do a quick validation
                    # In some cases, there really is a "negative" intergenic distance for overlapping CDS sequences
                    genic_distance = current_gene_start - upstream_gene_end
                    if genic_distance < 0 and current_gene_start < upstream_gene_start:
                        try:
                            ig_locus = '{0}--{1}_ig'.format(upstream_gene, locus)
                            genic_distance = agg_data[ig_locus]['Length']
                        except KeyError:
                            print("Caught something")
                            genic_distance = current_gene_start + (int(genome_length[chrom]['TotalLength']) - upstream_gene_end)


                # If the absolute upstream gene is on the positive strand, minus ig hops are antisense

                if upstream_gene_strand == 1:

                    #print("Upstream antisense for entry {0} - {1} on strand {2}".format(locus, upstream_gene, upstream_gene_strand))

                    antisense_data[locus+'_upstream'] = {'Reference': current_gene_reference,
                                                       'Locus': locus,
                                                         'Gene': current_gene_name,
                                                         'Start': current_gene_start,
                                                         'End': current_gene_end,
                                                         'Length': current_gene_length,
                                                       'AntisenseSites': current_gene_minus_sites,
                                                         'TotalSites': current_gene_total_sites,
                                                       'AntisenseCount': current_gene_minus_count,
                                                       'percentAntisenseCount': current_gene_percent_minus_count,
                                                       'TotalCount': current_gene_total_count,
                                                       'Antisense Contrib to Dval': current_gene_minus_dval,
                                                       'TotalDval': current_gene_dval,
                                                         'Product': current_gene_product,
                                                       'Downstream Locus': upstream_gene,
                                                       'Downstream Gene': upstream_gene_name,
                                                       'Downstream Locus Strand': upstream_gene_strand,
                                                         'Downstream Start': upstream_gene_start,
                                                         'Downstream End': upstream_gene_end,
                                                         'Downstream Length': upstream_gene_length,
                                                       'Intergenic Distance': genic_distance,
                                                       'Downstream PlusCount': upstream_gene_plus_count,
                                                       'Downstream MinusCount': upstream_gene_minus_count,
                                                       'Downstream DvalGenome': upstream_gene_dval,
                                                       'Downstream Product': upstream_gene_product
                                                       }

                # If the absolute upstream gene is on the negative strand, minus ig hops are sense
                if upstream_gene_strand == -1:

                    #print("Upstream sense for entry {0} - {1} on strand {2}".format(locus, upstream_gene, upstream_gene_strand))

                    sense_data[locus+'_upstream'] = {'Reference': current_gene_reference,
                                                           'Locus': locus,
                                                     'Gene': current_gene_name,
                                                     'Start': current_gene_start,
                                                     'End': current_gene_end,
                                                     'Length': current_gene_length,
                                                           'SenseSites': current_gene_minus_sites,
                                                     'TotalSites': current_gene_total_sites,
                                                     'SenseCount': current_gene_minus_count,
                                                           'percentSenseCount': current_gene_percent_minus_count,
                                                           'TotalCount': current_gene_total_count,
                                                           'Sense Contrib to Dval': current_gene_minus_dval,
                                                           'TotalDval': current_gene_dval,
                                                     'Product': current_gene_product,
                                                           'Downstream Locus': upstream_gene,
                                                           'Downstream Gene': upstream_gene_name,
                                                           'Downstream Locus Strand': upstream_gene_strand,
                                                     'Downstream Start': upstream_gene_start,
                                                     'Downstream End': upstream_gene_end,
                                                     'Downstream Length': upstream_gene_length,
                                                           'Intergenic Distance': genic_distance,
                                                           'Downstream PlusCount': upstream_gene_plus_count,
                                                           'Downstream MinusCount': upstream_gene_minus_count,
                                                           'Downstream DvalGenome': upstream_gene_dval,
                                                           'Downstream Product': upstream_gene_product
                                                           }

            if downstream_gene and 'End of chromosome' not in downstream_gene:

                downstream_gene_reference = agg_data[downstream_gene]['Reference']
                downstream_gene_name = agg_data[downstream_gene]['Gene']
                downstream_gene_strand = agg_data[downstream_gene]['Strand']
                downstream_gene_plus_sites = agg_data[downstream_gene]['PlusSites']
                downstream_gene_minus_sites = agg_data[downstream_gene]['MinusSites']
                downstream_gene_plus_count = agg_data[downstream_gene]['PlusCount']
                downstream_gene_minus_count = agg_data[downstream_gene]['MinusCount']
                downstream_gene_start = agg_data[downstream_gene]['Start']
                downstream_gene_end = agg_data[downstream_gene]['End']
                downstream_gene_length = agg_data[downstream_gene]['Length']
                downstream_gene_dval = agg_data[downstream_gene]['DvalGenome']
                downstream_gene_product = agg_data[downstream_gene]['Product']

                # intergenic distance hack
                if '_ig' in locus:
                    genic_distance = agg_data[locus]['Length']
                else:
                    # In case we are wrapping around the chromosome, do a quick validation
                    # In some cases, there really is a "negative" intergenic distance for overlapping CDS sequences
                    genic_distance = downstream_gene_start - current_gene_end
                    if genic_distance < 0 and downstream_gene_start < current_gene_start:
                        try:
                            ig_locus = '{0}--{1}_ig'.format(locus, downstream_gene)
                            genic_distance = agg_data[ig_locus]['Length']
                        except KeyError:
                            print("Caught something")
                            genic_distance = current_gene_start + (
                                    int(genome_length[chrom]['TotalLength']) - downstream_gene_end)

                # If the current gene and absolute downstream gene are sense, then it is a sense contribution
                if downstream_gene_strand == 1:

                    #print("Downstream sense for entry {0} - {1} on strand {2}".format(locus, downstream_gene, downstream_gene_strand))


                    sense_data[locus+'_downstream'] = {'Reference': current_gene_reference,
                                                       'Locus': locus,
                                                       'Gene': current_gene_name,
                                                       'Start': current_gene_start,
                                                       'End': current_gene_end,
                                                       'Length': current_gene_length,
                                                       'SenseSites': current_gene_plus_sites,
                                                       'TotalSites': current_gene_total_sites,
                                                       'SenseCount': current_gene_plus_count,
                                                       'percentSenseCount': current_gene_percent_plus_count,
                                                       'TotalCount': current_gene_total_count,
                                                       'Sense Contrib to Dval': current_gene_plus_dval,
                                                       'TotalDval': current_gene_dval,
                                                       'Product': current_gene_product,
                                                       'Downstream Locus': downstream_gene,
                                                       'Downstream Gene': downstream_gene_name,
                                                       'Downstream Locus Strand': downstream_gene_strand,
                                                       'Downstream Start': downstream_gene_start,
                                                       'Downstream End': downstream_gene_end,
                                                       'Downstream Length': downstream_gene_length,
                                                       'Intergenic Distance': genic_distance,
                                                       'Downstream PlusCount': downstream_gene_plus_count,
                                                       'Downstream MinusCount': downstream_gene_minus_count,
                                                       'Downstream DvalGenome': downstream_gene_dval,
                                                       'Downstream Product': downstream_gene_product
                                                       }

                # If the absolute downstream gene is on the negative strand, plus ig hops are antisense
                if downstream_gene_strand == -1:

                    #print("Downstream antisense for entry {0} - {1} on strand {2}".format(locus, downstream_gene, downstream_gene_strand))

                    antisense_data[locus+'_downstream'] = {'Reference': current_gene_reference,
                                                           'Locus': locus,
                                                           'Gene': current_gene_name,
                                                           'Start': current_gene_start,
                                                           'End': current_gene_end,
                                                           'Length': current_gene_length,
                                                           'AntisenseSites': current_gene_plus_sites,
                                                           'TotalSites': current_gene_total_sites,
                                                           'AntisenseCount': current_gene_plus_count,
                                                           'percentAntisenseCount': current_gene_percent_plus_count,
                                                           'TotalCount': current_gene_total_count,
                                                           'Antisense Contrib to Dval': current_gene_plus_dval,
                                                           'TotalDval': current_gene_dval,
                                                           'Product': current_gene_product,
                                                           'Downstream Locus': downstream_gene,
                                                           'Downstream Gene': downstream_gene_name,
                                                           'Downstream Locus Strand': downstream_gene_strand,
                                                           'Downstream Start': downstream_gene_start,
                                                           'Downstream End': downstream_gene_end,
                                                           'Downstream Length': downstream_gene_length,
                                                           'Intergenic Distance': genic_distance,
                                                           'Downstream PlusCount': downstream_gene_plus_count,
                                                           'Downstream MinusCount': downstream_gene_minus_count,
                                                           'Downstream DvalGenome': downstream_gene_dval,
                                                           'Downstream Product': downstream_gene_product
                                                           }

    return sense_data, antisense_data


# Function writes interaction file
def write_combo(outfile, payload, direction):

    fieldnames = []

    if direction == 'sense':
        fieldnames = ['Reference', 'Locus', 'Gene', 'SenseSites', 'TotalSites', 'SenseCount', 'percentSenseCount', 'TotalCount', 'Start', 'End', 'Length', 'Sense Contrib to Dval',
                      'TotalDval', 'Product', 'Downstream Locus', 'Downstream Gene', 'Downstream Locus Strand', 'Downstream Start', 'Downstream End', 'Downstream Length', 'Intergenic Distance',
                      'Downstream PlusCount', 'Downstream MinusCount', 'Downstream DvalGenome', 'Downstream Product']
    if direction == 'antisense':
        fieldnames = ['Reference', 'Locus', 'Gene', 'AntisenseSites', 'TotalSites', 'AntisenseCount', 'percentAntisenseCount', 'TotalCount', 'Start', 'End', 'Length', 'Antisense Contrib to Dval',
                      'TotalDval', 'Product', 'Downstream Locus', 'Downstream Gene', 'Downstream Locus Strand', 'Downstream Start', 'Downstream End', 'Downstream Length', 'Intergenic Distance',
                      'Downstream PlusCount', 'Downstream MinusCount', 'Downstream DvalGenome', 'Downstream Product']

    with open(outfile, 'w') as handle:

        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()

        for entry in sorted(payload):
            writer.writerow(payload[entry])

    return


# Function writes an aggregate file
def write_aggregate(outfile, payload):

    fieldnames = ['Reference', 'Locus', 'PlusSites', 'MinusSites', 'Sites', 'PlusCount', 'MinusCount', 'percentPlus', 'TotalCount', 'LargestHopContrib', 'Start', 'End',
                  'Strand',
                  'Length', 'ExpectedCount', 'DvalGenome', 'Gene', 'ProteinID', 'Product', 'Note']

    total_entries = 0
    total_cds = 0
    inserted_entry = 0
    total_inserted_cds = 0

    with open(outfile, 'w') as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter='\t')

        writer.writeheader()

        for entry in sorted(payload):

            total_entries += 1
            if int(payload[entry]['TotalCount']) > 0:
                inserted_entry += 1
            if '_ig' not in payload[entry]['Locus']:
                total_cds += 1
                if int(payload[entry]['TotalCount']) > 0:
                    total_inserted_cds += 1

            writer.writerow(payload[entry])

        agg_stat_payload = ["\nAGGREGATE STATS:",
                            "Total entries: {0}".format(total_entries),
                            "Total entries with insertion: {0} (Approximate total complexity: {1:.2f}%)".format(inserted_entry, 100*(inserted_entry/total_entries)),
                            "Total CDS: {0}".format(total_cds),
                            "Total CDS with insertion: {0} (Approximate CDS complexity: {1:.2f}%)\n".format(total_inserted_cds, 100*(total_inserted_cds/total_cds))]

        for stat in agg_stat_payload:
            print(stat)

    return agg_stat_payload


# Loads FASTQ files from dir
def load_fastq_files(fastq_dir):

    fastqs = []
    try:
        for fastq_file in os.listdir(fastq_dir):
            if 'fastq' in fastq_file.split('.'):
                fastqs.append(os.path.join(fastq_dir, fastq_file))
    except FileNotFoundError:
        print("FASTQ folder not detected! Please make a folder called 'fastq' in the hopcount directory and place your read files in that folder.")
        quit(2)

    return fastqs


# Finds a reference file
def find_ref():

    gbk_files = []
    for k in os.listdir(os.getcwd()):
        if k[-3:] == '.gb' or k[-4:] == '.gbk':
            print("Detected GenBank reference file: {0}".format(k))
            gbk_files.append(os.path.join(os.getcwd(), k))
    if len(gbk_files) > 1:
        print("Error! More than one reference GenBank file was detected in the hopcount directory. Please place only one file for mapping here.\n"
              "You may combine references in a combined GenBank file, such as two vibrio chromosomes, if you wish.")
        quit(4)
    if len(gbk_files) == 0:
        print("Error! No refernce Genbank file detected.")
        quit(4)
    else:
        return gbk_files[0]


# Loads runtime configuration
def load_hopcount_conf_file(conf_file):

    bowtie2_conf, hopcount_conf = dict(), dict()

    try:
        with open(conf_file, 'r') as f:
            for line in f:

                # this is a bowtie option
                if line[0] == '-':
                    bowtie2_conf[line.split(' ')[0]] = line.split(' ')[1][:-1]

                # this is a hopcount option
                if line[0:3] == 'var':
                    hopcount_conf[line.split('=')[0]] = line.split('=')[1][:-1]
    except FileNotFoundError:
        print("Warning! Configuration file expected at {0} not found. Please place the configuration file in the same directory as the hopcount script.".format(conf_file))
        quit()

    return bowtie2_conf, hopcount_conf


# makes a wiggle file from the transposition histogram
# wiggle format found here: https://useast.ensembl.org/info/website/upload/wig.html
# IGV does not currently support multiple track wiggle files - make separate ones for plus, minus, total
def make_wiggle_files(hist, outdir, base):

    wig_plus = os.path.join(outdir, base + '_plus.wig')
    wig_minus = os.path.join(outdir, base + '_minus.wig')
    wig_total = os.path.join(outdir, base + '_total.wig')

    with open(wig_plus, 'w') as plus:
        with open(wig_minus, 'w') as minus:
            with open(wig_total, 'w') as total:

                # write track def files
                plus.write('track type=wiggle_0 name="{0} hopcount plusCount" visibility=full color=255,0,0\n'.format(base))
                minus.write('track type=wiggle_0 name="{0} hopcount minusCount" visibility=full color=0,0,255\n'.format(base))
                total.write('track type=wiggle_0 name="{0} hopcount totalCount" visibility=full color=0,0,0\n'.format(base))

                for chrom in hist:
                    plus.write('variableStep chrom={0}\n'.format(chrom))
                    minus.write('variableStep chrom={0}\n'.format(chrom))
                    total.write('variableStep chrom={0}\n'.format(chrom))

                    for pos in hist[chrom]:
                        plusCount, minusCount = hist[chrom][pos][0], hist[chrom][pos][1]
                        totalCount = plusCount + minusCount

                        plus.write('{0} {1}\n'.format(pos, plusCount))
                        minus.write('{0} {1}\n'.format(pos, -int(minusCount)))
                        total.write('{0} {1}\n'.format(pos, totalCount))

    return


# writes an output configuration file to a folder
def write_config(config_file, payload):

    with open(config_file, 'w') as g:

        g.write('### HOPCOUNT.PY RUNTIME CONFIGURATION AND STATISTICS ###\n')

        for line in payload:
            g.write(line + '\n')
        g.write('\n### END OF FILE ###')
    return


def main(ref, reads, output, bt2ops, hopcops):

    """
    This is the main flow script of HOPCOUNT. This is what occurs in order:

    1. Extract CDS information
        a. Extract CDS information from input GBK file
        b. Create a list of all potential intergenic regions based on mined CDS information

    2. Map the sequencing reads to the reference genome using Bowtie2.
        a. Convert GBK file to FASTA
        b. Map sequencing reads to FASTA file using bowtie2 call, with specified bt2 options in hopcount.conf

    3. From the resulting SAM file, create a transposition histogram with CDS annotation
        a. For each mapped read in the SAM file, determine the position of the read, and bin into the appropriate region
            based on provided positional thresholding parameters
        b. Output the transposition histogram to file
        c. Output WIGGLE files to file (for genome browsing and CIRCOS diagrams)

    4. Using the transposition histogram, create an aggregate file in which individual entries in the histogram are
        aggregated based on CDS
        a. Aggregate reads based on read thresholding parameters provided
        b. Calculate various statistics, including complexity, DvalGenome, read contributions, etc.
        c. Output aggregate table to file

    5. Piece together synergistic and antagonistic interactions that are relevant when using inducible promoters that
        read off the transposon

    6. Finally, write runtime files and clean up
        a. Write runtime file, including basic stats, files, options, etc.
        b. Deal with the bulky SAM file as specified
        c. Remove temporary files generated, such as bt2 index

    This flow is iterated over every read file from steps 2b-6, since the genome is the same each time.

    :param ref: GBK file reference of the organism
    :param reads: Path to the reads
    :param output: Output directory
    :param bt2ops: Dictionary of bowtie2 options
    :param hopcops: Dictionary of hopcount runtime options
    :return: Nothing!
    """

    # 1a Load up the gbk file - extract annotations
    print("Extracting gene information from {0}...".format(ref))
    gene_info, genome_lengths, split_cds, cds_lengths = extract_gene_info(gbk=ref, is_circular=hopcops['var_circular_genome'])

    # Generate a thresholded gene list
    threshold_cds = threshold_gene_list(gene_info, percent_threshold=float(hopcops['var_percent_cutoff']),
                                        genome_lengths=genome_lengths, is_circular=hopcops['var_circular_genome'])

    # 1b Create a list of all possible intergenic entries based on CDS data
    print("Creating intergenic entries based on CDS data...")
    gene_info, ig_lengths = generate_ig_regions(gene_info, is_circular=hopcops['var_circular_genome'], genome_lengths=genome_lengths)

    all_lengths = {}
    for chromosome in cds_lengths:
        all_lengths[chromosome] = {**cds_lengths[chromosome], **ig_lengths[chromosome]}

    # 2a Convert gbk to fasta for bowtie2
    print("Converting {0} to FASTA file for bowtie2...".format(ref))
    ref_basename = os.path.splitext(os.path.basename(ref))[0]
    gbk_to_fasta(ref, ref_basename+'.fasta')

    # 2a Create index file
    print("Building index for mapping...")
    bowtie2_build(ref_basename+'.fasta', 'temp')

    # Iterate over read files
    for read_file in reads:

        runtime_config_payload = ['Read file: {0}'.format(read_file),
                                  'Reference file: {0}'.format(ref),
                                  'Bowtie2 configuration: {0}'.format(' '.join([item for k in bt2ops for item in (k, bt2ops[k])])),
                                  'Hopcount.py configuration: {0}\n'.format(' '.join([item for k in hopcops for item in (k, hopcops[k])]))]

        print("Processing read file {0}...".format(read_file))

        # Make output directory
        file_output_name = os.path.basename(read_file).split('.')[0]
        file_output_dir = os.path.join(output, file_output_name)
        if not os.path.exists(file_output_dir):
            os.mkdir(os.path.join(file_output_dir))
            print("Created output folder directory at {0}/".format(os.path.join(output, file_output_dir)))

        sam_file_loc = os.path.join(file_output_dir, file_output_name+'.sam')

        # 2b map reads to reference
        print("Mapping reads from {0} to {1}...".format(read_file, ref))
        bowtie2_stats = bowtie2_unpaired('temp', read_file, sam_file_loc, bt2ops)

        runtime_config_payload.append('Bowtie2 read statistics:\n{0}'.format(bowtie2_stats))

        if not bowtie2_stats:
            print("Bowtie2 failed! Quitting...")
            quit(5)

        # 3a Obtain a position histogram from the sam file
        print("Obtaining transposition read position histogram from sam file...")
        transposition_read_histogram, hopcount_stats = get_transp_read_pos(sam_file_loc, int(hopcops['var_quality']))

        runtime_config_payload.append(hopcount_stats)

        # Now we want output a sorted file, called output_hopcount.csv, that contains this histogram data
        # We can perform percent within cds thresholding here using our thresholded gene list - makes aggregate dead simple
        # 3b
        hopcount_file_path = os.path.join(file_output_dir, file_output_name+'_hopcount.tsv')
        print("Obtaining genomic context and writing hopcount data to {0}...".format(hopcount_file_path))
        write_stats = write_hopcount_file(hopcount_file_path, threshold_cds, transposition_read_histogram,
                                          gene_info, percent_cutoff=float(hopcops['var_percent_cutoff']),
                                          genome_lengths=genome_lengths)

        runtime_config_payload += write_stats

        # 3c Make WIGGLE files
        try:
            if hopcops['var_make_wiggle_file'] == 'true':
                print("Creating WIGGLE files...")
                make_wiggle_files(transposition_read_histogram, file_output_dir, file_output_name)
        except ValueError:
            pass

        # Now let's read that sorted file and do fitness calculations.
        # 4a
        agg_file_path = os.path.join(file_output_dir, file_output_name+'_aggregate.tsv')
        print("Reading and compiling data from {0}...".format(hopcount_file_path))
        aggregate_data, agg_stats = aggregate_hopcount(hopcount_file_path, genome_lengths, read_count_threshold=int(hopcops['var_hop_threshold']), gene_info=gene_info)
        runtime_config_payload.append(agg_stats)

        print("Aligning additional gene data to aggregate data...")
        align_gene_info(aggregate_data, all_lengths)

        print("Performing expectation analysis...")
        expect_val(aggregate_data, genome_lengths)

        # Process synergistic and antagonistic intergenic ptac interactions
        print("Analyzing effects of outward reading promoter of intergenic insertions...")
        syn, ant = agg_analyze_combos(aggregate_data, gene_info, genome_length=genome_lengths, is_circular=hopcops['var_circular_genome'])

        # Write the payload
        print("Outputting aggregate file to {0}...".format(agg_file_path))
        agg_stat_payload = write_aggregate(agg_file_path, aggregate_data)
        runtime_config_payload += agg_stat_payload

        # Write payload
        syn_file = os.path.join(file_output_dir, file_output_name+'_sense.tsv')
        ant_file = os.path.join(file_output_dir, file_output_name+'_antisense.tsv')
        print("Outputting data to {0} and {1}...".format(syn_file, ant_file))
        write_combo(syn_file, syn, direction='sense')
        write_combo(ant_file, ant, direction='antisense')

        # Write config file
        run_file_path = os.path.join(file_output_dir, file_output_name+'_runconfig.txt')
        print("Writing configuration stats file {0}...".format(run_file_path))
        write_config(run_file_path, runtime_config_payload)

        # Remove SAM, or convert BAM files
        if hopcops['var_keep_sam_file'] == 'no':
            print("Deleting SAM file...")
            os.remove(sam_file_loc)
        if hopcops['var_keep_sam_file'] == 'bam':
            print('Converting to BAM file...')
            bamify(sam_file_loc, os.path.join(file_output_dir, file_output_name+'.bam'))

    # Remove temp files
    print("Removing temporary index and fasta files...")
    files_in_dir = os.listdir(os.getcwd())
    for file in files_in_dir:
        if file.startswith('temp'):
            os.remove(os.path.join(os.getcwd(), file))
    return


if __name__ == '__main__':

    # VARS FOR FILE LOCATIONS - SHOULD ADAPT TO CWD
    hopcount_conf_file = os.path.join(os.getcwd(), 'hopcount.conf')     # Configuration file for hopcount
    fastq_folder = os.path.join(os.getcwd(), 'fastq')                   # FASTQ file folder location
    output_folder = os.path.join(os.getcwd(), 'output')                 # Base output folder

    # Read FASTQ
    print("Loading FASTQ files...", end='')
    fastq_files = load_fastq_files(fastq_folder)
    if not fastq_files:
        print("Error! No FASTQ files detected in {0}. Please place at least one read file in that folder.".format(fastq_folder))
        quit(3)
    else:
        print("loaded {0} FASTQ read file(s).".format(len(fastq_files)))
        for f in fastq_files:
            print(f)

    # Read Ref. Scan the dir for a gb or gbk file.
    print("Finding reference file...", end='')
    ref_file = find_ref()
    print("detected at {0}".format(ref_file))

    # Read Options
    print("Reading configuration file {0}...".format(hopcount_conf_file), end='')
    bowtie2_ops, hopcount_ops = load_hopcount_conf_file(hopcount_conf_file)
    if (not bowtie2_ops) or (not hopcount_ops):
        print("Error! A configuration file was found, but no options were loaded. Please check the conf file at {0}.".format(hopcount_conf_file))
        quit(3)
    else:
        print("loaded options.")

    print('Bowtie2 options: {0}'.format(' '.join([item for k in bowtie2_ops for item in (k, bowtie2_ops[k])])))
    print('Hopcount options: {0}'.format(' '.join([item for k in hopcount_ops for item in (k, hopcount_ops[k])])))

    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        print("Created output folder at {0}/".format(output_folder))
        os.mkdir(output_folder)

    main(ref=ref_file, reads=fastq_files, output=output_folder, bt2ops=bowtie2_ops, hopcops=hopcount_ops)

    print("Done!")
