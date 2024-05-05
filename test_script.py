import os
import pytest


from bio_files_processor import FastaRecord, OpenFasta
from ultimate_bioinf_tool import (AminoAcidSequence, DNASequence, RNASequence)


@pytest.fixture
def input_data_fasta():
    fasta_path = './example_data/example_fasta.fasta'
    return fasta_path


class TestDataValue:
    @pytest.mark.parametrize(
        ('value', 'target'),
        [
            (('', '', ''), False),
            (('GTD323452', '', ''), True),
            (('', '', 'ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG'), False)
        ],
    )

    def test_fastarecord_booling(self, value, target):
        result = FastaRecord(*value)
        assert bool(result) == target

    def test_openfasta_value(self, input_data_fasta):
        with OpenFasta(input_data_fasta) as fasta:
            record = fasta.read_record()
        
        with open(input_data_fasta, 'r') as file:
            line = file.readline().strip().strip('>')
            id_, description = line.split(' ', 1)
            sequence = file.readline().strip() + file.readline().strip()
            result = FastaRecord(id_, description, sequence)
    
        assert record == result

class TestDataTypes:
    @pytest.mark.parametrize(
        ('seq_type','value'),
        [
            (RNASequence, 'ATGC'),
            (DNASequence, 'AUGC')
        ],
    )

    def test_bsw_typecheck(self, seq_type, value):
        with pytest.raises(TypeError):
            seq = seq_type(value)

    @pytest.mark.parametrize(
        ('seq_type','value'),
        [
            (RNASequence, 'AUGC'),
            (DNASequence, 'ATGC')
        ],
    )

    def test_bsw_return_type(self, seq_type, value):
        complement_type = type(seq_type(value).complement())
        assert complement_type == type(seq_type(value))

    def test_openfasta_return_type(self, input_data_fasta):
        with OpenFasta(input_data_fasta) as fasta:
            record = fasta.read_record()
        assert isinstance(record, FastaRecord)


class TestProtocolImplementation:

    def test_openfasta_is_context_manager(self, input_data_fasta):
        fasta = OpenFasta(input_data_fasta)
        assert hasattr(fasta, '__enter__') and hasattr(fasta, '__exit__')
    
    def test_openfasta_io_on_exit(self, input_data_fasta):
        with pytest.raises(ValueError):
            with OpenFasta(input_data_fasta) as fasta:
                line = fasta.read_record()
            fasta.read_record()

    def test_openfasta_is_iterator(self, input_data_fasta):
        fasta = OpenFasta(input_data_fasta)
        assert hasattr(fasta, '__iter__') and hasattr(fasta, '__next__')
