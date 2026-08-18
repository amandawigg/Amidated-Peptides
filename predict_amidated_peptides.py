"""
Prediction of candidate C-terminally amidated peptides from a secretome.

Wiggenhorn et al., "Global profiling of peptide amidation to discover
bioactive fragments of the secretome."

Classically secreted proteins are scanned from the signal-peptide cleavage site
to the C-terminus for a site of PAM-mediated amidation. Two cases are searched:
a glycine-dibasic tripeptide within the sequence (GKR, GKK, GRK, GRR), and a
terminal GR or GK where the protein sequence ends.

For each site an N-terminal cleavage position is taken as the most 3' dibasic
pair (KR, KK, RK, RR) downstream of the signal peptide, or the annotated
signal-peptide cleavage site. The predicted peptide is the sequence between that
position and the residue preceding the glycine.

Input is a UniProt TSV export with columns:
    0 Entry   2 Gene names   4 Modified residue   7 Signal peptide
    8 Peptide   9 Sequence   10 Propeptide   11 Chain

Run from the repository root.
"""

import csv
import os

os.makedirs('results', exist_ok=True)

# ── Files ────────────────────────────────────────────────────────────────────
secretome = open('data/mouse_secretome.tsv', 'r')
output    = open('results/mouse_amidated_peptides.csv', 'w', newline='')

# ── Parameters ───────────────────────────────────────────────────────────────
AMIDATION_MOTIFS = ('GRR', 'GKR', 'GKK', 'GRK')  # glycine-dibasic tripeptides
DIBASIC_MOTIFS   = ('RR', 'KR', 'KK', 'RK')      # convertase cleavage pairs
MIN_LENGTH       = 3    # shortest peptide reported
MAX_LENGTH       = 15   # longest peptide reported (exclusive); raise for the
                        # full predicted set at any length
FLANK            = 5    # residues of flanking sequence reported either side

# csv.writer quotes fields containing commas, which UniProt annotation columns
# (peptide / propeptide / chain / modified residue) frequently contain.
writer = csv.writer(output)
writer.writerow([
    'peptide name', 'peptide sequence', 'uniprot accession', 'GN',
    'peptide length', 'start position', 'stop position',
    'around N-term', 'around C-term', 'signal start', 'signal peptide end',
    'peptide info', 'propeptide info', 'chain info', 'mod info',
])

# ── Read the secretome ───────────────────────────────────────────────────────
secreted_genes = {}   # accession -> [name, sequence, signal_end, peptide_info,
                      #               propeptide_info, chain_info, mod_info]

for line in secretome:
    line = line.strip().split('\t')
    if len(line) > 1 and line[0] != "Entry":
        name = line[2]
        if name == '':
            name = line[0]
        accession = line[0]
        sequence  = line[9]
        ptm       = line[7].replace('"', '')

        # Only proteins with an annotated signal peptide are considered.
        if ptm[:6] == 'SIGNAL':
            signal       = line[7].split(';')[0]
            signal_end   = str(signal.split('..')[1])   # kept as text: may be '?'
            peptide_info = line[8]
            propeptide_info = line[10] if len(line) > 10 else ''
            chain_info      = line[11] if len(line) > 11 else ''
            mod_info        = line[4]

            secreted_genes[accession] = [name, sequence, signal_end, peptide_info,
                                         propeptide_info, chain_info, mod_info]

# ── Predict peptides ─────────────────────────────────────────────────────────
pep_list = []   # sequences already written, so each is reported once

for gene in secreted_genes:
    name            = secreted_genes[gene][0]
    seq             = secreted_genes[gene][1]
    signal_end_aa   = secreted_genes[gene][2]
    peptide_info    = secreted_genes[gene][3]
    propeptide_info = secreted_genes[gene][4]
    chain_info      = secreted_genes[gene][5]
    mod_info        = secreted_genes[gene][6]

    count = 1   # peptide number within this protein, in N-to-C order

    for position in range(len(seq) - 1):

        # Parentheses group the motif test so that the position test applies
        # to all four motifs (and binds tighter than or).
        if (seq[position:position + 3] in AMIDATION_MOTIFS) and (position > 1):

            fragment = seq[:position]

            if signal_end_aa.isdigit():
                signal_end = int(signal_end_aa)
            else:
                signal_end = -1

            # Most 3' dibasic pair anywhere upstream of the amidation motif.
            indices = []
            for motif in DIBASIC_MOTIFS:
                indices.append(fragment.rfind(motif))
            index = max(indices)

            if index > signal_end:
                start = index + 2          # cleave immediately after the pair
            elif signal_end != -1:
                start = signal_end
            else:
                start = 0

            pot_hormone = fragment[start:]

            # Flanking context; max(0, ...) guards the sequence start.
            if start > FLANK:
                upstream = seq[start - FLANK - 1:start + FLANK]
            else:
                upstream = seq[:start + FLANK]

            if len(seq[position + 3:]) > FLANK:
                downstream = seq[max(0, position - 4):position + FLANK + 2]
            else:
                downstream = seq[max(0, position - 4):]

            if MIN_LENGTH <= len(pot_hormone) < MAX_LENGTH:
                if pot_hormone not in pep_list:
                    pep_list.append(str(pot_hormone))

                    GN = name.split(' ')[0]

                    # Peptides are numbered by position in the precursor: the
                    # second site in a protein becomes PEP-GENE-2. The counter
                    # is per UniProt entry.
                    if count > 1:
                        pep_name = 'PEP-' + str(GN) + '-' + str(count)
                    else:
                        pep_name = 'PEP-' + str(GN)

                    if start == signal_end:
                        signal_start = 'yes'
                    else:
                        signal_start = 'no'

                    writer.writerow([
                        pep_name,
                        pot_hormone,
                        gene,
                        name,
                        len(pot_hormone),
                        start + 1,                    # reported 1-based
                        len(pot_hormone) + start,
                        upstream,
                        downstream,
                        signal_start,
                        signal_end_aa,
                        peptide_info,
                        propeptide_info,
                        chain_info,
                        mod_info,
                    ])
                    count += 1


    # ── Terminal GR / GK ─────────────────────────────────────────────────────
    # The second case: the protein sequence ends in GR or GK, so the glycine is
    # followed by a single basic residue at the C-terminus.
    if seq.endswith('GR') or seq.endswith('GK'):
        position = len(seq) - 2          # index of the glycine
        fragment = seq[:position]

        if signal_end_aa.isdigit():
            signal_end = int(signal_end_aa)
        else:
            signal_end = -1

        indices = []
        for motif in DIBASIC_MOTIFS:
            indices.append(fragment.rfind(motif))
        index = max(indices)

        if index > signal_end:
            start = index + 2
        elif signal_end != -1:
            start = signal_end
        else:
            start = 0

        pot_hormone = fragment[start:]

        if start > FLANK:
            upstream = seq[start - FLANK - 1:start + FLANK]
        else:
            upstream = seq[:start + FLANK]
        downstream = seq[max(0, position - 4):]

        if MIN_LENGTH <= len(pot_hormone) < MAX_LENGTH:
            if pot_hormone not in pep_list:
                pep_list.append(str(pot_hormone))

                GN = name.split(' ')[0]
                if count > 1:
                    pep_name = 'PEP-' + str(GN) + '-' + str(count)
                else:
                    pep_name = 'PEP-' + str(GN)

                if start == signal_end:
                    signal_start = 'yes'
                else:
                    signal_start = 'no'

                writer.writerow([
                    pep_name, pot_hormone, gene, name, len(pot_hormone),
                    start + 1, len(pot_hormone) + start, upstream, downstream,
                    signal_start, signal_end_aa, peptide_info, propeptide_info,
                    chain_info, mod_info,
                ])
                count += 1

output.close()
secretome.close()
print('Peptides written:', len(pep_list))
