"""
Title: DEVELOPMENT OF BIOINFORMATICS TOOL FOR THE DESIGN OF AMPLICONS IN VIRUS 
SEQUENCING
Author: Enrique Caban Centeno
Institution: Master in Bioinformatics, University of Valencia

This script is part of my Master's thesis in Bioinformatics. The code
provides a tool for primer generation and genome sequence analysis. It
uses the Biopython library to manage and manipulate aligned genomic sequences,
generate consensus sequences, and design primers with specific properties.
The tool includes functions to evaluate and optimize primers, ensuring they
meet criteria for melting temperature and sequence compatibility.
Additionally, it offers functionalities to visualize primer positions and
generate detailed reports in various formats. The tool processes aligned
genome sequences from FASTA files, performs primer design in batches, and
visualizes the results, making it a versatile resource for genomic research
and bioinformatics applications.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import time


class GenomeMarker:
    def __init__(self, amplicon_length, overlap_length, sequences=None):
        self.genome_records = sequences if sequences else []
        self.primers = []
        self.primers_batch2 = []
        self.features = []
        self.probable_nucleotides = {}
        self.amplicon_length = amplicon_length
        self.overlap_length = overlap_length

    def generate_consensus_or_return_sequence(self):
        """
        Generates a consensus sequence if there are multiple genome records,
        otherwise returns the sequence of the single genome record.

        :return: Consensus sequence or single genome record sequence as a
        string.
        """
        
        if len(self.genome_records) > 1:
            return self.generate_consensus()
        
        else:
            return str(self.genome_records[0].seq)

    def generate_consensus(self):
        """
        Generates a consensus sequence from multiple genome records.

        :return: Consensus sequence as a string.
        """
        
        consensus_sequence = ""
        for i in range(len(self.genome_records[0].seq)):
            nucleotide_counts = {"a": 0, "c": 0, "g": 0, "t": 0, "n": 0, "-": 0}
            for record in self.genome_records:
                nucleotide = record.seq[i].lower()
                if nucleotide in nucleotide_counts:
                    nucleotide_counts[nucleotide] += 1
            consensus_nucleotide = self.calculate_consensus_nucleotide(
                nucleotide_counts, len(self.genome_records)
            )
            consensus_sequence += consensus_nucleotide
        return consensus_sequence
        
    def calculate_consensus_nucleotide(self, counts, total):
        """
        Calculates the consensus nucleotide based on the counts of each
        nucleotide in the given position across all genome records.

        :param counts: Dictionary with counts of each nucleotide.
        :param total: Total number of genome records.
        :return: Consensus nucleotide as a single character string.
        """
        iupac_ambiguity_codes = {
            frozenset(["a", "g"]): "r",
            frozenset(["c", "t"]): "y",
            frozenset(["g", "c"]): "s",
            frozenset(["a", "t"]): "w",
            frozenset(["g", "t"]): "k",
            frozenset(["a", "c"]): "m",
            frozenset(["c", "g", "t"]): "b",
            frozenset(["a", "g", "t"]): "d",
            frozenset(["a", "c", "t"]): "h",
            frozenset(["a", "c", "g"]): "v",
            frozenset(["a", "c", "g", "t"]): "n",
        }
        nucleotides_above_90 = {
            nt for nt, perc in counts.items() if (perc / total) * 100 > 90
        }
        if nucleotides_above_90:
            return nucleotides_above_90.pop()
        nucleotides_above_10 = {
            nt for nt, perc in counts.items() if (perc / total) * 100 > 10
        }
        if nucleotides_above_10:
            return iupac_ambiguity_codes.get(frozenset(nucleotides_above_10), "n")
        return "n"

    def save_sequences_with_consensus_to_fasta(self, output_filename):
        """
        Saves the generated consensus sequence and original sequences to a 
        FASTA file.

        :param output_filename: Name of the output FASTA file.
        """
        print("Generating the Consensus Sequence")
        main_sequence = self.generate_consensus_or_return_sequence()
        main_sequence_record = SeqRecord(
            Seq(main_sequence),
            id="Consensus_Sequence",
            description="Generated or single consensus sequence",
        )
        original_sequences = (
            [record for record in self.genome_records]
            if len(self.genome_records) > 1
            else []
        )
        all_sequences = [main_sequence_record] + original_sequences

        with open(output_filename, "w") as output_file:
            SeqIO.write(all_sequences, output_file, "fasta")

    def process(self, min_tm, tm_difference_threshold, start_position):
        """
        Processes the genome records to generate primers, saving the details
        to a text file.

        :param min_tm: Minimum annealing temperature (Tm) for primers.
        :param tm_difference_threshold: Maximum allowed difference in Tm
        between forward and reverse primers.
        :param start_position: Starting position for primer generation.
        """
        print("Desining the first batch the primers across the Consensus Genome")
        sequence_to_process = self.generate_consensus_or_return_sequence()
        pos = start_position  # Starting position for the first primer

        # While there is enough sequence left to process another amplicon
        while pos + self.amplicon_length <= len(sequence_to_process):
            # Generar el primer forward al inicio del amplicón
            forward_primer_start = sequence_to_process[
                pos : pos + 20
            ]  # Suponiendo longitud de primer de 17
            reverse_complement_start = str(
                Seq(forward_primer_start).reverse_complement()
            )

            # Generate the forward primer at the beginning of the amplicon
            reverse_primer_end = sequence_to_process[
                pos + self.amplicon_length - 20 : pos + self.amplicon_length
            ]
            reverse_complement_end = str(
                Seq(reverse_primer_end).reverse_complement()
            )

            # Generate the reverse primer at the end of the amplicon
            start_pos_forward = pos + 1
            end_pos_forward = pos + 20
            start_pos_reverse = pos + self.amplicon_length - 20 + 1
            end_pos_reverse = pos + self.amplicon_length

            # Adjust positions for forward and reverse primers
            self.primers.append(
                (
                    forward_primer_start,
                    "forward",
                    start_pos_forward,
                    end_pos_forward,
                    reverse_complement_end,
                    "reverse",
                    start_pos_reverse,
                    end_pos_reverse,
                )
            )

            # Move the position for the next amplicon
            pos = end_pos_reverse + (
                self.amplicon_length - (2 * self.overlap_length)
            )

        # Check primer positions to ensure they meet specifications.
        self.check_primer_positions(self.primers)

        # Call to check_tm_difference
        self.check_tm_difference(tm_difference_threshold, self.primers)

        # After generating all primers, call to check_spaces
        self.check_spaces(self.generate_consensus_or_return_sequence(), self.primers)
        # Optimize primers to meet specifications
        self.evaluate_min_tm_primers(min_tm, self.primers)
        self.optimize_and_check_primers(
            min_tm,
            tm_difference_threshold,
            self.primers,
            20,
        )
        # After generating all primers and before final optimization
        consensus_sequence = self.generate_consensus_or_return_sequence()
        self.extend_primers(
            min_tm, self.primers, consensus_sequence, "Batch1"
        )

        # Call to calculate_reverse_primer_degeneration for each primer
        for primer in self.primers:
            seq_primer_forward = primer[0]
            seq_primer_reverse = primer[4]
            self.calculate_primer_degneration(seq_primer_forward)
            self.calculate_primer_degneration(seq_primer_reverse)


    def calculate_tm(self, sequence):
        """
        Calculates the annealing temperature (Tm) of a DNA sequence
        considering ambiguous nucleotides.
        Uses average values for ambiguous nucleotides.

        :param sequence: DNA sequence for which Tm will be calculated.
        :return: Annealing temperature (Tm) as a float.
        """
        # Tm values for non-ambiguous nucleotides
        tm_values = {"A": 2, "T": 2, "C": 4, "G": 4}

        # Average contributions for ambiguous nucleotides based on the
        # possible nucleotides they represent
        ambiguous_tm = {
            "R": (tm_values["A"] + tm_values["G"]) / 2,
            "Y": (tm_values["C"] + tm_values["T"]) / 2,
            "S": (tm_values["G"] + tm_values["C"]) / 2,
            "W": (tm_values["A"] + tm_values["T"]) / 2,
            "K": (tm_values["G"] + tm_values["T"]) / 2,
            "M": (tm_values["A"] + tm_values["C"]) / 2,
            "B": (tm_values["C"] + tm_values["G"] + tm_values["T"]) / 3,
            "D": (tm_values["A"] + tm_values["G"] + tm_values["T"]) / 3,
            "H": (tm_values["A"] + tm_values["C"] + tm_values["T"]) / 3,
            "V": (tm_values["A"] + tm_values["C"] + tm_values["G"]) / 3,
            "N": (tm_values["A"] + tm_values["C"] + tm_values["G"] + tm_values["T"])
            / 4,
        }

        # Update the Tm values dictionary to include ambiguous nucleotides
        tm_values.update(ambiguous_tm)

        # Calculate the Tm of the sequence
        Tm = sum(tm_values[nt.upper()] for nt in sequence if nt.upper() in tm_values)

        return float(Tm)

    def check_primer_positions(self, primers):
        """
        Checks if the last 5 bases of forward primers and the first 5 bases of
        reverse primers are valid (i.e., contain only A, T, C, G).

        :param primers: List of tuples containing primer sequences and their
        start and end positions.
        :return: Boolean indicating whether all primer positions are valid.
        """
        all_positions_valid = True
        for i, (
            primer_forward,
            _,
            forward_start,
            forward_end,
            primer_reverse,
            _,
            reverse_start,
            reverse_end,
        ) in enumerate(primers):
            forward_valid = all(base in "ATCG" for base in primer_forward[-5:].upper())
            reverse_valid = all(base in "ATCG" for base in primer_reverse[:5].upper())
            if not forward_valid:
                all_positions_valid = False
            if not reverse_valid:
                all_positions_valid = False
        return all_positions_valid

    def check_tm_difference(self, threshold_difference, primers):
        """
        Checks if the temperature difference (Tm) between any two primers
        exceeds a specified threshold.

        :param threshold_difference: The maximum allowed difference in Tm
        between any two primers.
        :param primers: List of primer information tuples, each containing
        primer sequences.
        :return: Boolean indicating whether all Tm differences are within the
        specified threshold.
        """
        all_differences_valid = True
        for i, primer_info_i in enumerate(primers):
            seq_primer_i = primer_info_i[0]
            tm_i = self.calculate_tm(seq_primer_i)
            for j, primer_info_j in enumerate(primers):
                if i == j:
                    continue
                seq_primer_j = primer_info_j[0]
                tm_j = self.calculate_tm(seq_primer_j)
                difference = abs(tm_i - tm_j)
                if difference > threshold_difference:
                    all_differences_valid = False
        return all_differences_valid


    def check_spaces(self, sequence, primers):
        """
        Checks if any part of the primers overlaps with spaces (indicated by '-') in the given sequence.

        :param sequence: DNA sequence string that may contain spaces.
        :param primers: List of tuples containing primer sequences and their
        positions.
        :return: None, but sets a flag if a primer overlaps a space.
        """
        space = "-"
        for i, (
            seq_primer_forward,
            _,
            forward_start,
            forward_end,
            seq_primer_reverse,
            _,
            reverse_start,
            reverse_end,
        ) in enumerate(primers):
            # Check if any nucleotide of the forward primer falls into a space
            for pos in range(forward_start - 1, forward_end):
                if sequence[pos] == space:
                    # Indicate an overlap without printing
                    self.primer_overlaps_space = True
                    break

            # Check if any nucleotide of the reverse primer falls into a space
            for pos in range(reverse_start - 1, reverse_end):
                if sequence[pos] == space:
                    # Indicate an overlap without printing
                    self.primer_overlaps_space = True
                    break

    def evaluate_min_tm_primers(self, min_tm, primers):
        """
        Evaluates whether both the forward and reverse primers in each pair have a minimum melting temperature (Tm).

        :param min_tm: Minimum required Tm for each primer.
        :param primers: List of primer tuples containing sequences of forward and reverse primers.
        :return: None, but sets a flag if any primer does not meet the minimum Tm.
        """
        for i, primer in enumerate(primers):
            seq_primer_forward = primer[0]
            seq_primer_reverse = primer[4]
            tm_forward = self.calculate_tm(seq_primer_forward)
            tm_reverse = self.calculate_tm(seq_primer_reverse)

            if tm_forward < min_tm or tm_reverse < min_tm:
                # Indicate that a primer does not meet the minimum Tm without printing
                self.primer_below_min_tm = True

    def optimize_and_check_primers(
        self,
        min_tm,
        tm_difference_threshold,
        primers,
        movement_range=20,
    ):
        """
        Optimizes primers based on their Tm, position validity, and potential Tm differences.

        :param min_tm: Minimum acceptable Tm for primers.
        :param tm_difference_threshold: Maximum acceptable Tm difference between primers.
        :param primers: List of primers to be optimized.
        :param movement_range: Range in which primers can be shifted for optimization.
        :return: None, modifies the primers list.
        """
        print("Optimizing primers that do not meet specifications...")

        optimized_primers = []
        sequence_to_process = self.generate_consensus_or_return_sequence()

        def evaluate_primer_score(
            primer_seq,
            tm,
            tm_differences,
            position_validity,
            primer_type,
            start_pos,
            end_pos,
            pre_primer_sequence,
        ):
            score = 0
            # Temperature score calculation
            if MIN_TM <= tm <= MAX_TM:
                score += 15
            elif SECONDARY_MIN_TM <= tm < MIN_TM or MAX_TM < tm <= SECONDARY_MAX_TM:
                score += 0  # No bonus or penalty

            # TM differences score
            if all(diff <= tm_difference_threshold for diff in tm_differences):
                score += 10
            else:
                score -= 15

            # Position and dimer check scores
            if self.check_primer_positions([primer_info]):
                score += 100
            if not self.find_dimerizing_pairs([primer_info]):
                score += 40
            if not self.check_spaces(sequence_to_process, [primer_info]):
                score += 40

            # Ambiguity check in the primer
            ambiguous_count = sum(
                1 for nt in primer_seq if nt.upper() not in "ATCG"
            )
            score -= ambiguous_count * PENALTY_FACTOR

            # Ambiguity check in the pre-primer region
            ambiguous_count_pre = sum(
                1 for nt in pre_primer_sequence if nt.upper() not in "ATCG"
            )
            if ambiguous_count_pre > 0:
                penalty_pre = 50 * ambiguous_count_pre
                score -= penalty_pre

            # Validity check score
            if position_validity:
                score += 70

            return score

        MIN_TM = min_tm
        MAX_TM = min_tm + 5
        SECONDARY_MIN_TM = 50
        SECONDARY_MAX_TM = 70
        PENALTY_FACTOR = 100

        for primer_info in primers:
            (
                forward_primer,
                forward_type,
                start_pos_forward,
                end_pos_forward,
                reverse_primer,
                reverse_type,
                start_pos_reverse,
                end_pos_reverse,
            ) = primer_info

            best_score_forward = -float("inf")
            best_config_forward = None
            for displacement in range(-movement_range, movement_range + 1):
                new_start_pos_forward = max(1, start_pos_forward + displacement)
                new_end_pos_forward = (
                    new_start_pos_forward + len(forward_primer) - 1
                )
                pre_primer_sequence = sequence_to_process[
                    max(0, new_start_pos_forward - 6) : new_start_pos_forward - 1
                ]
                new_seq_primer_forward = sequence_to_process[
                    new_start_pos_forward - 1 : new_end_pos_forward
                ]
                new_tm_forward = self.calculate_tm(new_seq_primer_forward)
                tm_differences = [
                    abs(new_tm_forward - self.calculate_tm(other[0]))
                    for other in primers
                    if other != primer_info
                ]
                position_validity = all(
                    base in "ATCG" for base in new_seq_primer_forward[-5:].upper()
                )
                score_forward = evaluate_primer_score(
                    new_seq_primer_forward,
                    new_tm_forward,
                    tm_differences,
                    position_validity,
                    "Forward",
                    new_start_pos_forward,
                    new_end_pos_forward,
                    pre_primer_sequence,
                )

                if score_forward > best_score_forward:
                    best_score_forward = score_forward
                    best_config_forward = (
                        new_seq_primer_forward,
                        forward_type,
                        new_start_pos_forward,
                        new_end_pos_forward,
                    )

            best_score_reverse = -float("inf")
            best_config_reverse = None
            seq_primer_as_forward = str(Seq(reverse_primer).reverse_complement())
            start_pos_as_forward = (
                sequence_to_process.rfind(seq_primer_as_forward) + 1
            )
            if start_pos_as_forward == 0:
                raise ValueError(
                    "Converted forward primer not found in consensus sequence."
                )

            for displacement in range(-movement_range, movement_range + 1):
                new_start_pos_as_forward = max(
                    1, start_pos_as_forward + displacement
                )
                new_end_pos_as_forward = (
                    new_start_pos_as_forward + len(seq_primer_as_forward) - 1
                )
                pre_primer_sequence = sequence_to_process[
                    max(0, new_start_pos_as_forward - 6) : new_start_pos_as_forward - 1
                ]
                new_seq_primer_as_forward = sequence_to_process[
                    new_start_pos_as_forward - 1 : new_end_pos_as_forward
                ]
                new_tm_as_forward = self.calculate_tm(new_seq_primer_as_forward)
                tm_differences = [
                    abs(new_tm_as_forward - self.calculate_tm(other[4]))
                    for other in primers
                    if other != primer_info
                ]
                position_validity = all(
                    base in "ATCG" for base in new_seq_primer_as_forward[:5].upper()
                )
                score_reverse = evaluate_primer_score(
                    new_seq_primer_as_forward,
                    new_tm_as_forward,
                    tm_differences,
                    position_validity,
                    "Reverse",
                    new_start_pos_as_forward,
                    new_end_pos_as_forward,
                    pre_primer_sequence,
                )

                if score_reverse > best_score_reverse:
                    best_score_reverse = score_reverse
                    best_config_reverse = (
                        new_seq_primer_as_forward,
                        reverse_type,
                        new_start_pos_as_forward,
                        new_end_pos_as_forward,
                    )

            if best_config_forward and best_config_reverse:
                optimized_forward_seq = best_config_forward[0]
                optimized_forward_start = best_config_forward[2]
                optimized_forward_end = best_config_forward[3]
                reconverted_reverse_seq = str(
                    Seq(best_config_reverse[0]).reverse_complement()
                )
                reconverted_reverse_start = best_config_reverse[2]
                reconverted_reverse_end = best_config_reverse[3]

                optimized_primers.append(
                    (
                        optimized_forward_seq,
                        "forward",
                        optimized_forward_start,
                        optimized_forward_end,
                        reconverted_reverse_seq,
                        "reverse",
                        reconverted_reverse_start,
                        reconverted_reverse_end,
                    )
                )

        primers[:] = optimized_primers

    def extend_primers(self, min_tm, primers, consensus_sequence, file_prefix):
        """
        Extends primer sequences to reach or exceed a specified minimum Tm by
        adding nucleotides from a consensus sequence.

        :param min_tm: Minimum Tm that each primer should reach after
        extension.
        :param primers: List of primers to be extended.
        :param consensus_sequence: Consensus sequence from which nucleotides
        will be added to the primers.
        :param file_prefix: Prefix for the filenames where extension results
        will be logged.
        :return: None, modifies the primers list.
        """
        print("Extending primers to meet minimum Tm requirements...")

        for i, primer in enumerate(primers):
            (
                seq_primer_forward,
                forward_type,
                start_pos_forward,
                original_end_pos_forward,
                seq_primer_reverse,
                reverse_type,
                start_pos_reverse,
                end_pos_reverse,
            ) = primer

            tm_forward = self.calculate_tm(seq_primer_forward)
            tm_reverse = self.calculate_tm(seq_primer_reverse)
            seq_primer_forward_extended = seq_primer_forward
            seq_primer_reverse_extended = seq_primer_reverse

            # Extending forward primer
            while tm_forward < min_tm and start_pos_forward > 1:
                start_pos_forward -= 1
                nucleotide_to_add = consensus_sequence[start_pos_forward - 1]
                seq_primer_forward_extended = (
                    nucleotide_to_add + seq_primer_forward_extended
                )
                tm_forward = self.calculate_tm(seq_primer_forward_extended)

            # Convert reverse primer to forward equivalent to find correct positions in the consensus
            seq_primer_as_forward = str(
                Seq(seq_primer_reverse).reverse_complement()
            )
            # Finding the converted sequence in the consensus, preserving original reverse primer positions
            start_pos_as_forward = (
                consensus_sequence.rfind(seq_primer_as_forward) + 1
            )
            if start_pos_as_forward == 0:
                raise ValueError(
                    f"Converted forward primer not found in consensus sequence for primer {i+1}."
                )
            tm_as_forward = self.calculate_tm(seq_primer_as_forward)

            # Extend as a forward primer
            while tm_as_forward < min_tm and start_pos_as_forward > 1:
                start_pos_as_forward -= 1
                nucleotide_to_add = consensus_sequence[start_pos_as_forward - 1]
                seq_primer_as_forward = nucleotide_to_add + seq_primer_as_forward
                tm_as_forward = self.calculate_tm(seq_primer_as_forward)

            # Convert back to reverse after extension
            seq_primer_reverse_extended = str(
                Seq(seq_primer_as_forward).reverse_complement()
            )
            # Adjust start and end positions according to the original positions
            start_pos_reverse_extended = start_pos_reverse - (
                len(seq_primer_reverse_extended) - len(seq_primer_reverse)
            )
            end_pos_reverse_extended = (
                start_pos_reverse_extended + len(seq_primer_reverse_extended) - 1
            )
            tm_reverse = self.calculate_tm(
                seq_primer_reverse_extended
            )  # Recalculate the Tm for the extended reverse primer

            # Update the primer tuple in the list
            primers[i] = (
                seq_primer_forward_extended,
                forward_type,
                start_pos_forward,
                original_end_pos_forward,
                seq_primer_reverse_extended,
                reverse_type,
                start_pos_reverse_extended,
                end_pos_reverse_extended,
            )


    def check_primer_dimerization(self, primer1, primer2, min_overlap=10):
        """
        Checks if two primers can dimerize, considering only overlaps that meet or exceed a minimum length.

        :param primer1: Sequence of the first primer.
        :param primer2: Sequence of the second primer.
        :param min_overlap: Minimum overlap length to consider for
        dimerization.
        :return: Boolean indicating if the primers can dimerize.
        """
        primer1 = primer1.upper()
        primer2_complement = Seq(primer2.upper()).reverse_complement()

        if len(primer1) < min_overlap or len(primer2_complement) < min_overlap:
            return False

        for j in range(min_overlap, len(primer1) + 1):
            if primer1[-j:] == primer2_complement[:j]:
                return True
        return False

    def find_dimerizing_pairs(self, primers):
        """
        Finds all pairs of primers from a list that can potentially form
        dimers.

        :param primers: List of primers to check for potential dimerization.
        :return: List of tuples, each containing two primers that can
        potentially dimerize.
        """
        dimerizing_pairs = []
        for i, primer1_data in enumerate(primers):
            for j, primer2_data in enumerate(primers):
                if i >= j:  # Avoid duplicates and self-comparisons
                    continue
                if self.check_primer_dimerization(primer1_data[0], primer2_data[0]):
                    dimerizing_pairs.append((primer1_data, primer2_data))
        return dimerizing_pairs

    def calculate_primer_degneration(self, seq_primer):
        """
        Calculates the degeneration score of a primer based on predefined
        degeneration values for each type of nucleotide.

        :param seq_primer: The primer sequence for which to calculate
        degeneration.
        :return: The average degeneration score as a float.
        """

        inverse_degeneration_map = {
            "A": 4,
            "T": 4,
            "C": 4,
            "G": 4,
            "R": 3,
            "Y": 3,
            "S": 3,
            "W": 3,
            "K": 3,
            "M": 3,
            "B": 2,
            "D": 2,
            "H": 2,
            "V": 2,
            "N": 1,
        }
        # Calculate the degree of degeneration
        total_degeneration = 0
        for nucleotide in seq_primer.upper():
            total_degeneration += inverse_degeneration_map.get(
                nucleotide, 4
            )  # Assume non-degenerate as default value

        # Calculate the average degeneration
        average_degeneration = total_degeneration / len(seq_primer)

        

        return average_degeneration

    def process_batch_2(self, min_tm, tm_difference_threshold, start_position):
        """
        Processes the second batch of primers, starting from a specified
        position and applying various checks and optimizations.

        :param min_tm: Minimum melting temperature required for the primers.
        :param tm_difference_threshold: Maximum allowable difference in
        melting temperatures between primers.
        :param start_position: Position in the sequence from which to start
        generating primers.
        :return: None, but modifies internal state.
        """
        print("Processing the second batch of primers...")

        sequence_to_process = self.generate_consensus_or_return_sequence()
        pos = (
            start_position + self.amplicon_length - self.overlap_length
        )  # Adjust initial position for the second batch

        while pos + self.amplicon_length <= len(sequence_to_process):
            # Generate the forward primer at the beginning of the amplicon
            forward_primer_start = sequence_to_process[
                pos : pos + 20
            ]
            reverse_complement_start = str(
                Seq(forward_primer_start).reverse_complement()
            )

            # Generate the reverse primer at the end of the amplicon
            reverse_primer_end = sequence_to_process[
                pos + self.amplicon_length - 20 : pos + self.amplicon_length
            ]
            reverse_complement_end = str(
                Seq(reverse_primer_end).reverse_complement()
            )

            # Adjust positions for forward and reverse primers
            start_pos_forward = pos + 1
            end_pos_forward = pos + 20
            start_pos_reverse = pos + self.amplicon_length - 20 + 1
            end_pos_reverse = pos + self.amplicon_length

            # Save the primers in the primers list for the second batch
            self.primers_batch2.append(
                (
                    forward_primer_start,
                    "forward",
                    start_pos_forward,
                    end_pos_forward,
                    reverse_complement_end,
                    "reverse",
                    start_pos_reverse,
                    end_pos_reverse,
                )
            )

            # Move the position for the next amplicon
            pos = end_pos_reverse + (
                self.amplicon_length - (2 * self.overlap_length)
            )

        # Call to check primer positions for the second batch
        self.check_primer_positions(self.primers_batch2)

        # Call to check Tm difference for the second batch
        self.check_tm_difference(tm_difference_threshold, self.primers_batch2)
        self.check_spaces(
            self.generate_consensus_or_return_sequence(), self.primers_batch2
        )
        self.evaluate_min_tm_primers(min_tm, self.primers_batch2)
        self.optimize_and_check_primers(
            min_tm,
            tm_difference_threshold,
            self.primers_batch2,
            20,
        )
        consensus_sequence = (
            self.generate_consensus_or_return_sequence()
        )  # This line should already be in your method
        self.extend_primers(
            min_tm, self.primers_batch2, consensus_sequence, "Batch2"
        )  # Make sure to pass consensus sequence here too
        dimerizing_pairs_batch2 = self.find_dimerizing_pairs(self.primers_batch2)
        if dimerizing_pairs_batch2:
            print(
                f"Found {len(dimerizing_pairs_batch2)} pairs of primers that could dimerize."
            )
            for pair in dimerizing_pairs_batch2:
                print(f"Potential dimerization between {pair[0][0]} and {pair[1][0]}")
        else:
            print("No dimerizing pairs found.")


    def save_to_file(self, primers, filename, batch_description, standard_length=20):
        """
        Saves information about primers to a FASTA file, including detailed
        descriptions of each primer.

        :param primers: List of primer information to save.
        :param filename: The filename to save the FASTA records to.
        :param batch_description: Description of the batch for annotation
        purposes.
        :param standard_length: Standard length of primers for determining
        extensions.
        :return: None, outputs to a file.
        """
        fasta_records = []

        for i, (
            seq_primer_forward,
            forward_type,
            start_pos_forward,
            end_pos_forward,
            seq_primer_reverse,
            reverse_type,
            start_pos_reverse,
            end_pos_reverse,
        ) in enumerate(primers, start=1):
            Tm_forward = self.calculate_tm(seq_primer_forward)
            Tm_reverse = self.calculate_tm(seq_primer_reverse)

            # Calculate the effective extension based on the difference in length from the standard
            forward_extension_length = max(0, len(seq_primer_forward) - standard_length)
            reverse_extension_length = max(0, len(seq_primer_reverse) - standard_length)

            # Detailed descriptions with updated extension information
            forward_description = f"{batch_description} Forward Primer {i}; Pos: {start_pos_forward}-{end_pos_forward}; Tm: {Tm_forward:.2f}C"
            reverse_description = f"{batch_description} Reverse Primer {i}; Pos: {start_pos_reverse}-{end_pos_reverse}; Tm: {Tm_reverse:.2f}C"

            if forward_extension_length > 0:
                forward_description += (
                    f"; Ext: {forward_extension_length} bases at start"
                )
            if reverse_extension_length > 0:
                reverse_description += f"; Ext: {reverse_extension_length} bases at end"

            # Create and add FASTA records
            forward_record = SeqRecord(
                Seq(seq_primer_forward),
                id=f"F_Primer_{i}_{batch_description}",
                description=forward_description,
            )
            reverse_record = SeqRecord(
                Seq(seq_primer_reverse),
                id=f"R_Primer_{i}_{batch_description}",
                description=reverse_description,
            )

            fasta_records.extend([forward_record, reverse_record])

        # Save to FASTA file
        with open(filename, "w") as output_file:
            SeqIO.write(fasta_records, output_file, "fasta")

    def save_degeneration_scores_to_file(
        self, output_filename, primer_batch, batch_name
    ):
        """
        Saves the degeneration scores of primers to a file.

        :param output_filename: Filename where the degeneration scores will be saved.
        :param primer_batch: List of primers for which degeneration scores are calculated.
        :param batch_name: Name of the batch for annotation in the output file.
        :return: None, outputs to a file.
        """
        # 'w' to overwrite or create a new file
        with open(output_filename, "w") as output_file:
            output_file.write(f"Degeneration Scores for {batch_name}:\n")
            for primer in primer_batch:
                seq_primer_forward = primer[0]
                seq_primer_reverse = primer[4]
                forward_score = self.calculate_primer_degneration(seq_primer_forward)
                reverse_score = self.calculate_primer_degneration(seq_primer_reverse)
                output_file.write(
                    f"Forward Sequence: {seq_primer_forward}, Degeneration Score: {forward_score:.2f}\n"
                )
                output_file.write(
                    f"Reverse Sequence: {seq_primer_reverse}, Degeneration Score: {reverse_score:.2f}\n"
                )

    def calculate_amplicons(self):
        """
        Calculates the positions of amplicons based on primer information from both batches.

        :return: A sorted list of tuples, each representing the start and end positions of an amplicon.
        """
        amplicons = []
        # Calculate amplicons for the first batch
        for i in range(0, len(self.primers)):
            start_amplicon = self.primers[i][2]  # start position of the forward primer
            end_amplicon = self.primers[i][7]  # end position of the reverse primer
            amplicons.append((start_amplicon, end_amplicon))

        # Calculate amplicons for the second batch similarly
        for i in range(0, len(self.primers_batch2)):
            start_amplicon = self.primers_batch2[i][2]
            end_amplicon = self.primers_batch2[i][7]
            amplicons.append((start_amplicon, end_amplicon))

        # Sort the amplicons by their start position
        amplicons.sort(key=lambda x: x[0])
        return amplicons

    def generate_genome_diagram(self, start_position):
        """
        Generates a genome diagram indicating the positions of primers and amplicons.

        :param start_position: The starting position in the genome for the diagram.
        :return: None, outputs a graphical representation of primers and amplicons.
        """
        gd_diagram = GenomeDiagram.Diagram("Primer Scheme")
        gd_track_primers = gd_diagram.new_track(
            1, name="Primers", greytrack=False, scale=0.5
        )
        gd_feature_set = gd_track_primers.new_set()

        # Calculate amplicons and add them to the diagram
        amplicons = self.calculate_amplicons()
        for start, end in amplicons:
            adjusted_start = max(start, start_position)
            adjusted_end = max(end, start_position)
            if adjusted_start > adjusted_end:
                adjusted_start, adjusted_end = adjusted_end, adjusted_start

            # Use daltonismo-friendly colors
            if any(
                primer[2] == start and primer[1] == "forward" for primer in self.primers
            ):
                strand = +1
                color = colors.orange
            elif any(
                primer[2] == start and primer[1] == "forward"
                for primer in self.primers_batch2
            ):
                strand = -1
                color = colors.skyblue
            else:
                strand = None  # Or some default value

            if strand:
                feature = SeqFeature(
                    FeatureLocation(adjusted_start, adjusted_end, strand), strand=strand
                )
                gd_feature_set.add_feature(
                    feature, color=color, label=False, label_size=5, label_angle=0
                )

        # Adjust the primer colors for visibility
        for primer in self.primers:
            forward_feature = SeqFeature(
                FeatureLocation(primer[2], primer[3]), strand=+1
            )
            gd_feature_set.add_feature(
                forward_feature,
                name=f"{primer[2]}-{primer[3]}",
                label=True,
                color=colors.darkblue,
                label_position="start",
                label_size=6,
                label_angle=0,
                sigil="BOX",
                arrowshaft_height=1,
            )
            reverse_feature = SeqFeature(
                FeatureLocation(primer[6], primer[7]), strand=+1
            )
            gd_feature_set.add_feature(
                reverse_feature,
                name=f"{primer[6]}-{primer[7]}",
                label=True,
                color=colors.darkgreen,
                label_position="end",
                label_size=6,
                label_angle=0,
                sigil="BOX",
                arrowshaft_height=1,
            )

        for primer in self.primers_batch2:
            forward_feature = SeqFeature(
                FeatureLocation(primer[2], primer[3]), strand=-1
            )
            gd_feature_set.add_feature(
                forward_feature,
                name=f"{primer[2]}-{primer[3]}",
                label=True,
                color=colors.darkblue,
                label_position="start",
                label_size=6,
                label_angle=180,
                sigil="BOX",
                arrowshaft_height=1,
            )
            reverse_feature = SeqFeature(
                FeatureLocation(primer[6], primer[7]), strand=-1
            )
            gd_feature_set.add_feature(
                reverse_feature,
                name=f"{primer[6]}-{primer[7]}",
                label=True,
                color=colors.darkgreen,
                label_position="end",
                label_size=6,
                label_angle=180,
                sigil="BOX",
                arrowshaft_height=1,
            )

        # Draw and save the diagram
        gd_diagram.draw(
            format="linear", orientation="landscape", pagesize="A4", fragments=5
        )
        gd_diagram.write("primers_diagram.pdf", "PDF")

    def generate_degeneration_visualization_batch1(self):
        """
        Generates a visualization of degeneration scores for the first batch of primers.

        :return: None, outputs a graphical representation of degeneration scores.
        """
        # Visualization for the first batch of primers, including both forward and reverse primers
        degenerations = []
        labels = []  # Labels will now include 'F' or 'R' to indicate primer type
        for i, primer in enumerate(self.primers):
            forward_degeneration = self.calculate_primer_degneration(
                primer[0]
            )  # seq_primer_forward
            reverse_degeneration = self.calculate_primer_degneration(
                primer[4]
            )  # seq_primer_reverse
            degenerations.extend([forward_degeneration, reverse_degeneration])
            labels.extend(
                [f"F{i+1}", f"R{i+1}"]
            )  # We label forward primers with 'F' and reverse with 'R'
        self._generate_degeneration_graph(
            degenerations,
            labels,
            "Primer Degeneration Visualization - Batch 1",
            "degeneration_visualization_batch1.png",
        )

    def generate_degeneration_visualization_batch2(self):
        """
        Generates a visualization of degeneration scores for the second batch of primers.

        :return: None, outputs a graphical representation of degeneration scores.
        """
        # Visualization for the second batch of primers, including both forward and reverse primers
        degenerations = []
        labels = []  # Labels will now include 'F' or 'R' to indicate primer type
        for i, primer in enumerate(self.primers_batch2):
            forward_degeneration = self.calculate_primer_degneration(primer[0])
            reverse_degeneration = self.calculate_primer_degneration(primer[4])
            degenerations.extend([forward_degeneration, reverse_degeneration])
            labels.extend(
                [f"F{i+1}", f"R{i+1}"]
            )  # We label forward primers with 'F' and reverse with 'R'
        self._generate_degeneration_graph(
            degenerations,
            labels,
            "Primer Degeneration Visualization - Batch 2",
            "degeneration_visualization_batch2.png",
        )

    def _generate_degeneration_graph(self, degenerations, labels, title, filename):
        """
        Helper function to generate a bar graph of degeneration scores.

        :param degenerations: List of degeneration scores.
        :param labels: Labels corresponding to each score in 'degenerations'.
        :param title: Title for the graph.
        :param filename: Filename for saving the graph.
        :return: None, outputs a graphical representation to a file.
        """
        norm = plt.Normalize(min(degenerations), max(degenerations))
        cmap = plt.get_cmap("coolwarm_r")
        colors = cmap(norm(degenerations))
        plt.figure(figsize=(10, 4))
        plt.bar(range(len(degenerations)), [1] * len(degenerations), color=colors)
        plt.xticks(range(len(degenerations)), labels, rotation=45, ha="right")
        plt.yticks([])
        cbar = plt.colorbar(
            plt.cm.ScalarMappable(norm=norm, cmap=cmap), orientation="vertical"
        )
        cbar.set_label("Degeneration (4=Higher, 1=Lower)")
        plt.title(title)
        plt.tight_layout()
        plt.savefig(filename)
        plt.close()

    def generate_primers_excel(self):
        """
        Generates an Excel file containing detailed information about primers from all batches.

        :return: None, outputs an Excel file with primer details.
        """
        standard_length = 20  # Longitud estándar para comparar
        batches = [("Batch 1", self.primers), ("Batch 2", self.primers_batch2)]
        for batch_name, batch in batches:
            data = []  # Lista para almacenar la información de cada fila

            for primer in batch:
                (
                    seq_primer_forward,
                    _,
                    start_pos_forward,
                    end_pos_forward,
                    seq_primer_reverse,
                    _,
                    start_pos_reverse,
                    end_pos_reverse,
                ) = primer

                # Calcular la información de extensión para el primer forward
                extension_forward_info = "No extension"
                if len(seq_primer_forward) > standard_length:
                    extended_nucleotides_forward = seq_primer_forward[:3]
                    extension_positions_forward = (
                        f"{start_pos_forward}-{start_pos_forward+2}"
                    )
                    extension_forward_info = f"Extension: {extended_nucleotides_forward} (Positions: {extension_positions_forward})"

                # Calcular la información de extensión para el primer reverse
                extension_reverse_info = "No extension"
                if len(seq_primer_reverse) > standard_length:
                    extended_nucleotides_reverse = seq_primer_reverse[-3:]
                    extension_positions_reverse = (
                        f"{end_pos_reverse-2}-{end_pos_reverse}"
                    )
                    extension_reverse_info = f"Extension: {extended_nucleotides_reverse} (Positions: {extension_positions_reverse})"

                # Cálculo de Tm y puntuaciones de degeneración
                Tm_forward = self.calculate_tm(seq_primer_forward)
                Tm_reverse = self.calculate_tm(seq_primer_reverse)
                degeneration_forward = self.calculate_primer_degneration(
                    seq_primer_forward
                )
                degeneration_reverse = self.calculate_primer_degneration(
                    seq_primer_reverse
                )

                # Agregando la fila al conjunto de datos
                data.append(
                    [
                        seq_primer_forward,
                        f"{start_pos_forward}-{end_pos_forward}",
                        extension_forward_info,
                        Tm_forward,
                        degeneration_forward,
                        seq_primer_reverse,
                        f"{start_pos_reverse}-{end_pos_reverse}",
                        extension_reverse_info,
                        Tm_reverse,
                        degeneration_reverse,
                    ]
                )

            # Creación del DataFrame
            df = pd.DataFrame(
                data,
                columns=[
                    "Forward Primer",
                    "Forward Positions",
                    "Forward Extension",
                    "Forward Tm",
                    "Forward Degeneration Score",
                    "Reverse Primer",
                    "Reverse Positions",
                    "Reverse Extension",
                    "Reverse Tm",
                    "Reverse Degeneration Score",
                ],
            )

            # Guardar a Excel
            excel_filename = f"primers_{batch_name.replace(' ', '_').lower()}.xlsx"
            df.to_excel(excel_filename, index=False)


def main():
    # Ask the user for the desired FASTA file
    fasta_filename = input("Enter the name of the desired FASTA file: ")

    # Ask the user for the minimum desired Tm for the primers
    min_tm = float(
        input("Enter the minimum annealing temperature (Tm) desired for the primers: ")
    )

    # Ask the user for the maximum threshold of Tm difference allowed
    tm_difference_threshold = float(
        input("Enter the maximum threshold of Tm difference allowed between primers: ")
    )

    # Preguntar al usuario por la longitud del amplicón deseada
    amplicon_length = int(input("Enter the desired amplicon length: "))

    # Preguntar al usuario por la longitud del solapamiento deseado
    overlap_length = int(input("Enter the desired overlap length: "))

    start_position = int(input("Enter the start position in the genome: "))
    # Mark the start of execution time
    start_time = time.time()
    try:
        records = list(SeqIO.parse(fasta_filename, "fasta"))
        genome_marker = GenomeMarker(
            sequences=records,
            amplicon_length=amplicon_length,
            overlap_length=overlap_length,
        )
        genome_marker.save_sequences_with_consensus_to_fasta("output.fasta")
        genome_marker.process(min_tm, tm_difference_threshold, start_position)
        genome_marker.process_batch_2(min_tm, tm_difference_threshold, start_position)

        # Save the primers from the first batch in their own FASTA file
        genome_marker.save_to_file(
            genome_marker.primers, "primers_batch1.fasta", "Batch1"
        )
        # Save the primers from the second batch in another FASTA file
        genome_marker.save_to_file(
            genome_marker.primers_batch2, "primers_batch2.fasta", "Batch2"
        )

        # Save degeneration scores to separate files for each batch
        genome_marker.save_degeneration_scores_to_file(
            "degeneration_scores_batch1.txt", genome_marker.primers, "Batch 1"
        )
        genome_marker.save_degeneration_scores_to_file(
            "degeneration_scores_batch2.txt", genome_marker.primers_batch2, "Batch 2"
        )

        # Generate genome diagram
        genome_marker.generate_genome_diagram(start_position)
        # Generate degeneration visualizations for each batch
        genome_marker.generate_degeneration_visualization_batch1()
        genome_marker.generate_degeneration_visualization_batch2()
        # Generate Excel files for each batch of primers
        genome_marker.generate_primers_excel()

        print("Process completed successfully!")
    except FileNotFoundError:
        print("The specified file was not found.")

    # Mark the end of execution time
    end_time = time.time()

    # Calculate and display the total execution time
    total_time = end_time - start_time
    print(f"Total execution time: {total_time:.2f} seconds")


if __name__ == "__main__":
    main()
