

#secretome = open('mouse_secretome.tsv','r') ##Mouse
#output = open('Mouse_amidated_peptides.csv', 'w') ##Mouse

secretome = open('human_secretome.tsv', 'r') ##Human
output = open('Human_amidated_peptides.csv', 'w') ##Human



frag_dict={}

output_dict = {}
name_dict = {}


index = 0
secreted_genes = {}
gene_num = 1
for line in secretome:
    index +=1
    line = line.strip().split('\t')
    if len(line) > 1  and line[0] != "Entry":
        name = line[2]
        if name == '':
            name = line[0]
        accession = line[0]
        sequence = line[9]
        ptm = line[7].replace('"','')

        if ptm[:6] == 'SIGNAL':

            signal = line[7].split(';')[0]
            signal_end = str(signal.split('..')[1])
            peptide_info = line[8]
            if len(line)>10:
                propeptide_info = line[10]
            else:
                propeptide_info = ''
            if len(line)>11:
                chain_info = line[11]
            else:
                chain_info = ''
            mod_info = line[4]

            secreted_genes[accession] = [name, sequence, signal_end, peptide_info, propeptide_info, chain_info, mod_info]
            gene_num += 1


output.write('peptide name' +',' + 'peptide sequence' + ',' + 'uniprot accession' + ',' + 'GN' + ',' + 'peptide length' + ',' + 'start position' + ',' + 'stop position' + ',' + 'around N-term' + ',' + 'around C-term' + ',' + 'signal start' + ',' + 'signal peptide end' + ',' + 'peptide info' + ',' +  'propeptide info' + ',' + 'chain info' + ',' + 'mod info' + '\n')

pep_list = []

pep_num = 0
gene_list_pep = []
pep_gene_num = 0
count = 0
for gene in secreted_genes:
    seq = secreted_genes[gene][1]
    count = 1
    signal_end_aa = secreted_genes[gene][2]
    name = secreted_genes[gene][0]
    peptide_info = secreted_genes[gene][3]
    propeptide_info = secreted_genes[gene][4]
    chain_info = secreted_genes[gene][5]
    mod_info = secreted_genes[gene][6]

    for position in range(len(seq) - 1):
        if (seq[position:position+3] == 'GRR')or (seq[position:position+3] =='GKR') or (seq[position:position+3] =='GKK') or (seq[position:position+3] =='GRK') and (position > 1):


            fragment = seq[:position]

            if signal_end_aa != '?':
                signal_end = int(signal_end_aa)
            else:
                signal_end = -1
            indices = []
            start = 0
            dibasic = ['RR', 'KR', 'KK', 'RK']
            for motif in dibasic:
                pos_dibasic = fragment.rfind(motif, start)
                indices.append(pos_dibasic)
            index = max(indices)

            if index > signal_end:
                start = index +2
            elif signal_end != -1:
                start = signal_end
            else:
                start = 0
            pot_hormone = fragment[start:]
            if start > 5:

                upstream = seq[start - 6:start + 5]
            else:
                upstream = seq[:start + 5]
            if len(seq[position + 3:]) > 5:
                downstream = seq[position - 4:position + 7]
            else:
                downstream = seq[position - 4:]
            if len(pot_hormone) > 2 and len(pot_hormone) <31:
                if pot_hormone not in pep_list:
                    pep_list += [str(pot_hormone)]
                    pep_num += 1
                    num_pep = pep_num

                    if gene not in gene_list_pep:
                        gene_list_pep += [gene]
                    GN = name.split(' ')[0]

                    if num_pep < 10:
                        num_pep = '00' + str(num_pep)
                    elif num_pep < 100:
                        num_pep = '0' + str(num_pep)
                    if count > 1:
                        pep_name = 'PEP-' + str(GN) + ' ' + str(count)
                    else:
                        pep_name = 'PEP-' + str(GN)
         
                    if int(start) == int(signal_end_aa):
                        signal_start = 'yes'
                        
                    else:
                        signal_start = 'no'
           
                    output.write(str(pep_name) + ',' + str(pot_hormone) + ',' + gene + ',' + name + ',' + str(
                        len(pot_hormone)) + ',' + str(start + 1) + ',' + str(
                        len(pot_hormone) + start) + ',' + upstream + ',' + downstream + ',' + str(
                        signal_start) + ',' + str(
                        signal_end_aa) + ',' + peptide_info + ',' + propeptide_info + ',' + chain_info + ',' + mod_info + '\n')
                    count += 1




output.close()
secretome.close()

