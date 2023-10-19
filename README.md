# HW 6. Files
> *This is the repo for the sixth homework of the BI Python 2023 course*

This repository contains two scripts and dependent modules.

## `ultimate_bioinf_tool.py` script

### `run_dna_rna_tools` function

#### Title

`run_dna_rna_tools` function performs 4 basic transformation for nucleotide sequences. `run_dna_rna_tools` uses 6 secondary functions listed in `/modules/dna_rna_tools.py`.

##### Secondary functions

- `check_seq` - checks whether passed sequence is valid, i.e. is non-empty sequence that can be classified as DNA or RNA; if sequence is valid, returns positive exit code and corresponding complementary rule, otherwise returns negative exit code and `None`;
- `print_result` - prints brief information about success in processing sequences: number of successfuly processed sequences, number of non-valid sequences and pseudographical representation of their position in input.

All functions listed below take single sequence and complemetary rule for it and return processed sequence:

- `transcribe` — returns transcribed sequence, if RNA sequence is passed, returns sequence itself; takes complemetary rule to detect RNA sequences.
- `reverse` — returns reversed sequence; takes complemetary rule to provide call unification inside primary function (i.e. all functions were able to recieve 2 arguments).
- `complement` — returns complement sequence;
- `reverse_complement` — returns reversed complement sequence;

##### Primary function

`run_dna_rna_tools` takes arbitrary number of sequences terminated with desired operation and proceeds all given sequences one-by-one using secondary functions. If sequence is recognized as non-valid, appends sequence number with sequence itself to `corrupt_seqs` array. Returns tuple containing two list:

- `result` - list with operation results for each valid sequence;
- `corrupt_seqs` - list with non-valid sequences and their indices;


#### Installation

```
from ultimate_bioinf_tool import run_dna_rna_tools
```

#### Usage

Pass desired sequences terminated by desired operation. All argument must be passed separately, do not encapsulate them in array:

```
>>> run_dna_rna_tools('AtG', 'AUTG', 'GATTACA', 'aug', 'cguuc', '', 'transcribe')

Note: you have passed RNA sequence(s) for transcription. Output will contain RNA sequences themselves.

Processing result: [+-+++-]

4 sequence(s) out of 6 given have been processed successfully.
2 has been recognized as corrupted, i.e. non-DNA as well as non-RNA.
(['AuG', 'GAUUACA', 'aug', 'cguuc'], [(1, 'AUTG'), (5, '')])
```

### `run_protein_analyzer_tool` function

#### Title

`protein_analyzer_tool.py` function performs 5 operation on protein sequences. Sequence can be written both using 1 and 3-letter abbreviations. `protein_analyzer_tool.py` uses 10 secondary functions listed in `/modules/protein_analyzer_tool.py`.

##### Secondary functions

- `check_and_procees_seq` - checks whether the string is valid protein sequence, i.e. non-empty sequence composed of commonly accepted aminoacids abbreviation. Transforms valid sequence from 3 to 1-letter abbreviation if required
- `mann_whitney_u`- Mann-Whitney U-test. Used to compare aminoacids composition in given sequence with average protein composition provided by Uniprot. Used as a second step of checkup in `check_and_procees_seq` function if sequence uses 1-letter abbreviation.;
- `decomposition` - decomposes 3 letter-sequence into list of single aminoacids. Used in `check_and_procees_seq` function if sequence uses 3-letter abbreviation.
- `seq_transform` - transforms aminoacid abbreviation format from 3-letter into 1-letter. Used in `check_and_procees_seq` function if sequence uses 3-letter abbreviation.

- `aa_content_check` - returns aminoacids content of the protein.
- `seq_length` - returns length of the protein.
- `protein_mass` - returns molecular weight of the protein in g/mol.
- `protein_formula` - returns molecular formula of the protein
- `aa_chain_charge` - returns charge of the protein (pH=7).

- `print_result` - Prints brief information about success in processing sequences: number of successfuly processed sequences, number of non-valid sequences and pseudographical representation of their position in input.


##### Primary function

- `run_protein_analyzer_tool` - main function in module. Takes arbitrary number of sequences terminated with desired operation and proceeds all given sequences one-by-one using secondary functions. If sequence is recognized as non-valid, appends sequence number with sequence itself to `corrupt_seqs` array. Function returns array containing proceeded sequences (or sequence itself if result consists of single value). To get access to `corrupt_seqs` array uncomment line 81.

#### Installation

```
from ultimate_bioinf_tool import run_protein_analyzer_tool
```

#### Usage

Run `run_protein_analyzer_tool` function. This function provides interface for all 5 operations from `OPERATIONS` dictionary. Takes various number of positional arguments and one keyword-only argument:

- First `n` arguments - protein sequences;
- Latter positional argument - desired operation from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge";
- `abbrevition` keyword-only argument. Should be type integer, 1 for 1-letter abbreviation and 3 for 3-letter.

Returns tuple containing two list:

- `result` - list with operation results for each valid sequence;
- `corrupt_seqs` - list with non-valid sequences and their indices;


Get molecular mass for insulin:

```
>>> run_protein_analyzer_tool("AMALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN", "protein_mass", abbreviation=1)
All 1 sequence(s) processed successfully

(14090, [])
```


Get aminoacids content for various peptide:


```
>>> run_protein_analyzer_tool("AsnAspAspAsn", "content_check", abbreviation=3)
All 1 sequence(s) processed successfully

({'A': 0.0,
  'R': 0.0,
  'N': 50.0,
  'D': 50.0,
  'C': 0.0,
  'Q': 0.0,
  'E': 0.0,
  'G': 0.0,
  'H': 0.0,
  'I': 0.0,
  'L': 0.0,
  'K': 0.0,
  'M': 0.0,
  'F': 0.0,
  'P': 0.0,
  'S': 0.0,
  'T': 0.0,
  'W': 0.0,
  'Y': 0.0,
  'V': 0.0},
 [])
```


Get protein mass:


```
>>> run_protein_analyzer_tool("", "protein_mass", abbreviation=1)
Processing result: [-]

0 sequence(s) out of 1 given have been processed successfully.
1 has been recognized as corrupted, i.e. non-protein

([], (0, ''))
```


#### Troubleshooting

`run_protein_analyzer_tool` raises errors in two cases:

*   Operation is not one from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge". If you are sure that input is correct, perform spell check.
*   Argument for `abrreviation` parameter is not integer from 1 or 3.


In other cases `run_protein_analyzer_tool` will not halt the execution. Troubleshooting can be performed with second element in tuple returned by `run_protein_analyzer_tool`, `corrupt_seqs` list. This list contains sequences recognized as non-valid together with their indices in original sequence, `(<sequence_index>, <sequence>)`. Sequence is suggested to be non-valid in these cases:


*   If sequence is not type `str`. Other iterable objects are not supported by the time.
*   Sequence is empty string.
*   Sequence contain symbols except 20 commonly accepted aminoacid abbeviations.

### `run_fastq_processor.py` function

#### Title

`run_fastq_processor` filters reads presented in input fasta file into specified output file using three metrics:
-   GC-content;
-   sequence length;
-   average phred quality.
First function select reads into dictionary using metrics and then writes filtered reads to specified output file.

##### Secondary functions

- `process_paths` - сhecks input path and generates name and folder for output file. If no `output_filename` specified, takes name of input file. If no `output_path` specified, stores output file in `result` folder. Note that passing empty string as name for folder will make function to store output file in the same directory as input file; if output filename will be not specified either, input file will be overwriten.
- `process_file` - reads fasta file and transforms its content into dictionary in form of reads names as keys and tuple with reads array and phred quality array as values.
- `save_output` - writes elements dictionary containing reads names as keys and tuple with reads array and phred quality array as values into file.
- `is_in_range` - checks whether value belongs to desired interval. Used in `check_seq_and_bounds`.
- `check_seq_and_bounds` - сounts GC-content and average quality, if sequence length do not exceed length bounds - checks whether both of them belong to specified intervals. If sequence does, returns `True`, otherwise returns `False`.

##### Primary function

`run_fastq_processor` filters reads presented in input fasta file into dictionary using three metrics, GC-content,sequence length and average phred quality. Then writes filtered reads to file named as specified by `output_filename` and stores it in `output_path` folder inside directory with input file.

#### Installation

```
from ultimate_bioinf_tool import run_fastq_processor
```

#### Usage

If pass more than `input_path`, better use keyword arguments:

```
fq_input = 'D:/.../example_data/example_fastq.fastq'

run_fastq_processor(input_path = fq_input, output_filename = 'out_fastq')
```

## `bio_files_processor.py`

### Title

#### Secondary functions

- `write_seqs_dict` - consequently writes each dictionary value to file as individual line in file.
- `clear_seqs_dict` - clears dictionary values if values are lists.

#### Primary functions

- `convert_multiline_fasta_to_oneline` - converts multiline fasta sequence to oneline.
- `change_fasta_start_pos` - Shifts starting nucleotide position to desired position for one line fasta file.
- `select_genes_from_gbk_to_fasta` - in progress.

### Installation

```
import bio_files_processor
```

### Usage

```
fa_multiline_input = 'D:/.../example_data/example_multiline_fasta.fasta'

# >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
# ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
# GGTGCTGTG

convert_multiline_fasta_to_oneline(fa_multiline_input, 'out_oneline_fasta')

# >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
# ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG

```

```
>>> fa_input = 'D:/Users/.../example_data/example_oneline_fasta.fasta
# >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
# ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT

>>> change_fasta_start_pos(fa_input, 2, 'out_oneline_fasta')

# >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
# GGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTAC
```

## Contacts

We hope our module provides useful tool for your work. If you encounter any errors, please mail one from our team: 

Belikova Angelina - kiit@gmail.com
Implemented: `protein_mass`, `seq_length`.

Aryuna Ayusheeva - aryuna.ayusheeva.1998@mail.ru
Implemented: `aa_content_check`, `aa_chain_charge`.

Bredov Denis - d2707bredov@gmail.com
Teamlead. Implemented: `ultimate_bioinf_tool.py` and `bio_files_processor.py` scripts, `dna_rna_tools` and `fastq_processor` modules and `mann_whitney_u`, `decomposition`, `seq_transform`, `protein_formula` and `check_and_procees_seq` functions from `protein_analyzer_tool` module.
