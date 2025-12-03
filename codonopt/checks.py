import pandas as pd
from Bio.Seq import Seq
import subprocess
import tempfile

# Restriction sites
RESTRICTION_SITES = ['GGTCTC', 'GAGACC', 'CGTCTC', 'GAGACG']  # BsaI, BsmBI

# Canonical SD motifs
SD_PATTERNS = ["AGGAGG", "GGAGG", "AGGGA", "AGGAG", "GGAG", "AGGA", "GAGG"]


def check_translation(df, dna_col='DNA sequence', prot_col='Sequence'):
	"""Return a DataFrame of rows where the DNA does NOT match the protein."""
	errors = []
	for idx, row in df.iterrows():
		dna_seq = row[dna_col].upper()
		expected_prot = row[prot_col]

		# translate, using TAG for * (GCE)
		translated = str(Seq(dna_seq).translate(table=1, to_stop=False))

		# In GCE, all stops (*) are TAG codons; translation will produce * as expected
		# Compare sequences
		if translated != expected_prot:
			errors.append((idx, expected_prot, translated))
	return pd.DataFrame(errors, columns=['ID','Expected','Translated'])


def check_restriction_sites(df, dna_col='DNA sequence', sites=RESTRICTION_SITES):
	"""Return a DataFrame of rows that contain any restriction site."""
	hits = []
	for idx, row in df.iterrows():
		dna_seq = row[dna_col].upper()
		found_sites = [site for site in sites if site in dna_seq]
		if found_sites:
			hits.append((idx, found_sites))
	return pd.DataFrame(hits, columns=['ID','RestrictionSites'])


def check_sd_sites(df, dna_col='DNA sequence', sd_motifs=SD_PATTERNS):
	"""Return a DataFrame of rows that contain any canonical SD motif."""
	hits = []
	for idx, row in df.iterrows():
		dna_seq = row[dna_col].upper()
		found_sd = [sd for sd in sd_motifs if sd in dna_seq]
		if found_sd:
			hits.append((idx, found_sd))
	return pd.DataFrame(hits, columns=['ID','SDMotifs'])


def check_amber_stop_codons(df, dna_col='DNA sequence', prot_col='Sequence'):
	"""
	Check that all '*' (amber stop) in protein sequence correspond to TAG codons in DNA.
	Returns a DataFrame of any violations.
	"""
	# Only sequences with '*' need checking
	df_with_amber = df[df[prot_col].str.contains(r"\*")]

	hits = []

	for idx, row in df_with_amber.iterrows():
		dna_seq = row[dna_col].upper()
		protein_seq = row[prot_col]

		# Split DNA into codons
		codons = [dna_seq[i:i+3] for i in range(0, len(dna_seq), 3)]

		for pos, aa in enumerate(protein_seq):
			if aa == "*":
				codon = codons[pos]
				if codon != "TAG":
					hits.append((idx, pos, codon))

	return pd.DataFrame(hits, columns=['ID', 'ProteinPos', 'CodonFound'])


def check_functional_sd_sites(df, dna_col='DNA sequence', sd_motifs=SD_PATTERNS):
	"""
	Return a DataFrame of rows that contain an SD motif 3–15 bp upstream
	of a valid bacterial start codon (ATG, GTG, TTG).
	"""
	
	# Build regex: SD motif + spacer (3–15 bp) + start codon
	sd_regex = re.compile(
		rf"({'|'.join(sd_motifs)}).{{3,13}}({'|'.join(START_CODONS)})",
		re.IGNORECASE
	)

	hits = []

	for idx, row in df.iterrows():
		dna_seq = row[dna_col].upper()
		
		for m in sd_regex.finditer(dna_seq):
			sd_seq = m.group(1)
			start_codon = m.group(2)
			sd_pos = m.start(1)
			start_pos = m.start(2)
			spacer_len = start_pos - (sd_pos + len(sd_seq))

			hits.append({
				'ID': idx,
				'SDMotif': sd_seq,
				'StartCodon': start_codon,
				'SD_Position': sd_pos,
				'Start_Position': start_pos,
				'Spacer_bp': spacer_len
			})

	return pd.DataFrame(hits)

