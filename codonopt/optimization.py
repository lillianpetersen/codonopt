from dnachisel import *
from dnachisel import SequencePattern
from dnachisel import AvoidHairpins
from Bio.Seq import Seq
from .utils import (
	reverse_translate,
	get_stop_codon_constraints,
)

def codonopt_ecoli(protein_seq, stop_codon="TAG", avoid_patterns=None, verbose=False):
	# starting DNA
	initial_dna = reverse_translate(protein_seq)

	# Make AvoidPattern() objects
	if avoid_patterns is None:
		avoid_patterns = ["BsaI_site", "Esp3I_site", "7xC", "7xG"]
	avoid_constraints = [AvoidPattern(p) for p in avoid_patterns]

	# avoid Shine–Dalgarno (SD) sequences
	SD_PATTERNS = ["AGGAGG", "GGAGG", "AGGGA", "AGGAG", "GGAG", "AGGA", "GAGG"]
	sd_constraints = [AvoidPattern(p) for p in SD_PATTERNS]

	# Find positions of stop codons (*) and enforce TAG (or other specified codon)
	stop_constraints = get_stop_codon_constraints(
		protein_seq,
		stop_codon=stop_codon
	)

	problem = DnaOptimizationProblem(
		sequence=initial_dna,
		constraints=[
			*avoid_constraints,
			*sd_constraints,
			*stop_constraints, 
			EnforceGCContent(mini=0.30, maxi=0.65, window=50),
			EnforceTranslation(translation=protein_seq),
		],
		objectives=[
			CodonOptimize(species='e_coli'),
			AvoidHairpins(stem_size=8, hairpin_window=200),
		],
	) 
	
	try:
		problem.resolve_constraints()
		problem.optimize()

		if verbose:
			print(problem.constraints_text_summary())
			print(problem.objectives_text_summary())

	except NoSolutionError:
		print('Initial optimization failed. Relaxing AGGA constraint')

		# remove the AGGA constraint
		SD_PATTERNS = ["AGGAGG", "GGAGG", "AGGGA", "AGGAG", "AAGGA", "GGAG", "GAGG"]
		sd_constraints = [AvoidPattern(p) for p in SD_PATTERNS]

		# if AGGA must be used, make sure it isn't 4-13bp upstream of ATG or GTG
		sd_regex = "AGGA[ATGC]{4,13}(ATG|GTG)"
		sd_pattern = SequencePattern(sd_regex)

		problem = DnaOptimizationProblem(
			sequence=initial_dna,
			constraints=[
				*avoid_constraints,
				*sd_constraints,
				*stop_constraints, 
				AvoidPattern(sd_pattern),
				EnforceGCContent(mini=0.30, maxi=0.65, window=50),
				EnforceTranslation(translation=protein_seq),
			],
		objectives=[
			CodonOptimize(species='e_coli'),
			AvoidHairpins(stem_size=8, hairpin_window=200),
		],
		) 

		problem.resolve_constraints()
		problem.optimize()

		if verbose:
			print(problem.constraints_text_summary())
			print(problem.objectives_text_summary())

	
	return problem.sequence


def codonopt_human(protein_seq, avoid_patterns=None, verbose=False):
	initial_dna = reverse_translate(protein_seq)

	# Make AvoidPattern() objects
	if avoid_patterns is None:
		avoid_patterns = ["BsaI_site", "Esp3I_site", "7xC", "7xG"]
	avoid_constraints = [AvoidPattern(p) for p in avoid_patterns]

	problem = DnaOptimizationProblem(
		sequence=initial_dna,
		constraints=[
			*avoid_constraints,
			AvoidRareCodons(species='h_sapiens', min_frequency=0.1),
			EnforceGCContent(mini=0.25, maxi=0.60),
			EnforceGCContent(mini=0.20, maxi=0.75, window=50),
			EnforceTranslation(translation=protein_seq),
		],
		objectives=[
			CodonOptimize(species='h_sapiens', method='match_codon_usage'),
			AvoidHairpins(stem_size=8, hairpin_window=200),
		],
	) 
	
	problem.resolve_constraints()
	problem.optimize()

	if verbose:
		print(problem.constraints_text_summary())
		print(problem.objectives_text_summary())


	return problem.sequence
