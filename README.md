# codonopt

Simple python package for codon-optimizing protein sequences for E. coli or human expression using DNA Chisel.

I wrote this code so we can all avoid Shine–Dalgarno sequences when optimizing for E. coli expression and save ourselves from many unfortunate internal start sites :)

Degenerate Shine-Dalgarno sequences include any of the following motifs:
- AGGAGG, GGAGG, AGGGA, AGGAG, GGAG, AGGA, GAGG

Functionality includes:
- Reverse translation and codon optimization for E. coli and human
- For E. coli: avoids Shine-Dalgarno sequences
- Avoids unwanted sequence motifs (default is BsaI, BsmbI, homopolymers)
- Avoids strong internal hairpins
- Constrains stop to defined codon (default is TAG, for GCE applications)
- GC content windows and global GC constraints

### Installation

From local source:
```
bash pip install .
```

From GitHub:
```
pip install git+https://github.com/lillianpetersen/codonopt.git
```


### Usage

```
from codonopt import codonopt_ecoli, codonopt_human, checks

protein = "SQPQKGRKPRDLELPLSPSLLGGPGPERTPGSGSGSGLQAPGPALTPSLLPTHTLTPVLLTPSSLPPSIHFWSTLSPIAPRSPAKLSFQFPSSGSAQVHIPSISVDGLSTPVVLSPGPQKP"
dna = codonopt_ecoli(protein)

protein_GCE = "SQPQKGRKPRDLELPLSPSLLGGPGPER*PGSGSGSGLQAPGPAL*PSLLPTHTL*PVLLTPSSLPPSIHFWSTLSPIAPR*PAKLSFQFPSSGSAQVHIPSISVDGLS*PVVL*PGPQKP"
dna = codonopt_ecoli(protein_GCE, stop_codon='TAG')

library['DNA sequence'] = library['Sequence'].apply(codonopt_ecoli)

translation_errors = checks.check_translation(library, dna_col='DNA sequence')
restriction_hits = checks.check_restriction_sites(library, dna_col='DNA sequence')
sd_hits = checks.check_functional_sd_sites(library, dna_col='DNA sequence')
amber_violations = checks.check_amber_stop_codons(library, dna_col='DNA sequence')

etc
```

