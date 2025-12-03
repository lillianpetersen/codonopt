from Bio.Data import CodonTable
from dnachisel import EnforceSequence

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

# Basic codon lookup table
aa_to_codon = {aa: codon for codon, aa in standard_table.forward_table.items()}

def reverse_translate(protein_seq, stop_codon="TAG"):
	dna_seq = []
	for aa in protein_seq:
		if aa == "*":
			dna_seq.append(stop_codon)
		else:
			dna_seq.append(aa_to_codon[aa])
	return "".join(dna_seq)

def get_stop_codon_constraints(protein_seq, stop_codon="TAG"):
	constraints = []
	for i, aa in enumerate(protein_seq):
		if aa == "*":
			pos = i * 3
			constraints.append(
				EnforceSequence(sequence=stop_codon, location=(pos, pos + 3))
			)
	return constraints
