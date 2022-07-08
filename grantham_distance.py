import itertools
import os

import numpy as np
import pandas as pd

########################################################################################################################
### System Constants ###
########################################################################################################################


class HLAAlleleSequenceMapping:
    def __init__(self, HLA_fas_path: str):
        self.load_HLA_allele_AA_mapping(HLA_fas_path)

    def parse_fas_section(self, section: str):
        """"""
        assert ">" not in section
        allele, allele_sequence = section.split('\n', maxsplit=1)
        return f"{allele[0]}*{allele[1:3]}:{allele[3:]}", allele_sequence.replace('\n', '')

    def load_HLA_allele_AA_mapping(self, HLA_fas_path: str) -> str:
        """ """
        with open(HLA_fas_path, 'r') as file:
            self.HLA_allele_AA_mapping = dict(self.parse_fas_section(section) for section in file.read().split(">")[1:])

        return self.HLA_allele_AA_mapping

    def _map(self, alleles: str):
        if isinstance(alleles, str):
            return self.HLA_allele_AA_mapping.get(alleles, None)

        return [self._map(allele) for allele in alleles]


class GranthamDistance:
    def __init__(self, grantham_AA_table_path: str):
        self.AA_distance_table = pd.read_csv(grantham_AA_table_path, sep='\t', index_col=0)

        self.AA_distance_dict = {}
        for i in self.AA_distance_table.index:
            for j in self.AA_distance_table.columns:
                self.AA_distance_dict[i + j] = self.AA_distance_table.loc[i, j]

    def __repr__(self):
        return str(self.AA_distance_table)

    def AA_pair_distance(self, AA_1: str, AA_2: str) -> int:
        """ """
        return self.AA_distance_dict.get(AA_1.upper() + AA_2.upper(), None)

    def sequence_pair_distance(self, AA_sequence_1: str, AA_sequence_2: str) -> float:
        """"""
        if not (AA_sequence_1 and AA_sequence_2):
            return None

        assert len(AA_sequence_1) == len(AA_sequence_2), f""
        normalize_length = len(AA_sequence_1)

        distances = [self.AA_pair_distance(AA_1, AA_2) for AA_1, AA_2 in zip(AA_sequence_1, AA_sequence_2)
                     if AA_1 != AA_2]

        if any(dist is None for dist in distances):
            return None

        return np.sum(distances) / normalize_length

    def sequence_group_distance(self, sequences: list, mean: bool = True) -> float:
        """ Mean of all the pairwise distances"""

        sequence_pairs = itertools.combinations(sequences, 2)
        distances = [self.sequence_pair_distance(*pair) for pair in sequence_pairs]
        if any(dist is None for dist in distances):
            return None

        return np.mean(distances)


########################################################################################################################
### End ###
########################################################################################################################
