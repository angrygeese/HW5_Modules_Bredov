# HW 5. Modules
> *This is the repo for the fifth homework of the BI Python 2023 course*

# `dna_rna_tools.py` description

This module contains `run_dna_rna_tools` function that performs 4 basic transformation for nucleotide sequences. `run_dna_rna_tools` uses 6 secondary functions.

## Secondary functions

- `check_seq` - checks whether passed sequence is valid, i.e. can be classified as DNA or RNA; if sequence is valid, returns positive exit code and corresponding complementary rule, otherwise returns negative exit code and `None`;
- `print_result` - prints brief information about success in processing sequences: number of successfuly processed sequences, number of non-valid sequences and pseudographical representation of their position in input.

All functions listed below take single sequence and complemetary rule for it and return processed sequence:

- `transcribe` — returns transcribed sequence, if RNA sequence is passed, returns sequence itself; takes complemetary rule to detect RNA sequences.
- `reverse` — returns reversed sequence; takes complemetary rule to provide call unification inside primary function (i.e. all functions were able to recieve 2 arguments).
- `complement` — returns complement sequence;
- `reverse_complement` — returns reversed complement sequence;

## Primary function

- `run_dna_rna_tools` - main function in module. Takes arbitrary number of sequences terminated with desired operation and proceeds all given sequences one-by-one using secondary functions. If sequence is recognized as non-valid, appends sequence number with sequence itself to `corrupt_seqs` array. By default, function returns array containing proceeded sequences (or sequence itself if result consists of single value). To get access to `corrupt_seqs` array uncomment line 81.

### `dna_rna_tools.py` usage

Just pass desired sequences first and then operation to the function. All argument must be passed separately, do not encapsulate them in array.

# `protein_analyzer_tool.py`

## Title

This module contains the `protein_analyzer_tool.py` function that performs 5 operation on protein sequences. Operation are maintained using 5 secondary functions.


## Usage

To run the `protein_analyzer_tool.py`, first import it as module


```
import protein_analyzer_tool
```


and run its main function, `run_protein_analyzer_tool` function. This function provides interface for all 5 operations from `OPERATIONS` dictionary. Takes various number of positional arguments and one keyword-only argument:

- First `n` arguments - protein sequences;
- Latter positional argument - desired operation from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge";
- `abbrevition` keyword-only argument. Should be type integer, 1 for 1-letter abbreviation and 3 for 3-letter.

Returns tuple containing two list:

- `result` - list with operation results for each valid sequence;
- `corrupt_seqs` - list with non-valid sequences and their indices;


## Examples

Get molecular mass in g/mol for insulin:

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


Such case provides no result:


```
>>> run_protein_analyzer_tool("", "protein_mass", abbreviation=1)
Processing result: [-]

0 sequence(s) out of 1 given have been processed successfully.
1 has been recognized as corrupted, i.e. non-protein

([], (0, ''))
```
Explanation of the latter is provided in next section


## Troubleshooting

`run_protein_analyzer_tool` raises errors in two cases:

*   Operation is not one from list: "content_check", "seq_length", "protein_formula", "protein_mass", "charge". If you are sure that input is correct, perform spell check.
*   Argument for `abrreviation` parameter is not integer from 1 or 3.


In other cases `run_protein_analyzer_tool` will not halt the execution. In other scenarios troubleshooting can be performed using second element in tuple returned by `run_protein_analyzer_tool`, `corrupt_seqs` list. This list contains sequences recognized as non-valid together with their indices in original sequence. in form of tuple `(<sequence_index>, <sequence>)`. Sequence is suggested to be non-valid in these cases:


*   If sequence is not type `str`. Other iterable objects are not supported by the time.
*   Sequence is empty string.

# Contacts
We hope our module provides useful tool for your work. If you encounter any errors, please mail one from our team: 

Belikova Angelina - kiit@gmail.com
Implemented: `protein_formula`, `protein_mass`, `seq_length`.

Aryuna Ayusheeva - aryuna.ayusheeva.1998@mail.ru
Implemented: `aa_content_check`, `aa_chain_charge`.

Bredov Denis - d2707bredov@gmail.com
Teamlead. Implemented: `dna_rna_tools` module and `Mann_Whitney_U`, `decomposition`, `seq_transform`, `check_and_procees_seq`, `print_result`, `run_protein_analyzer_tool` from `` module.
