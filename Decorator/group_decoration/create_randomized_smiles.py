#!/usr/bin/env python

import argparse
import os

import utils.log as ul
import utils.chem as uc
import utils.scaffold as usc
import utils.spark as us


def parse_args():
    """Parses input arguments."""
    parser = argparse.ArgumentParser(
        description="Creates many sets with a given seed.")
    parser.add_argument("--input-smi-path", "-i",
                        help="Path to a SMILES file to convert with scaffolds and decorations.", type=str, required=True)
    parser.add_argument("--output-smi-folder-path", "-o",
                        help="Path to a folder that will have the converted SMILES files.", type=str, required=True)
    parser.add_argument("--num-files", "-n",
                        help="Number of SMILES files to create (numbered from 000 ...) [DEFAULT: 1]",
                        type=int, default=1)
    parser.add_argument("--num-partitions", "-p", help="Number of SPARK partitions to use [DEFAULT: 1000]",
                        type=int, default=1000)

    return parser.parse_args()


def _to_sliced_mol(row):
    scaffold_smi, decorations, _ = row.split("\t")
    decoration_smis = decorations.split("|")
    return usc.SlicedMol(uc.to_mol(scaffold_smi), {i: uc.to_mol(dec) for i, dec in enumerate(decoration_smis)})


def _format_training_set_row(sliced_mol):
    scaff_smi, dec_smis = sliced_mol.to_smiles(variant="random")
    return (scaff_smi, usc.ATTACHMENT_SEPARATOR_TOKEN.join(dec_smis))


def main():
    """Main function."""
    args = parse_args()

    sliced_mols_rdd = SC.textFile(args.input_smi_path) \
        .repartition(args.num_partitions) \
        .map(_to_sliced_mol)\
        .persist()

    os.makedirs(args.output_smi_folder_path, exist_ok=True)

    format_func = _format_training_set_row

    for i in range(args.num_files):
        with open("{}/{:03d}.smi".format(args.output_smi_folder_path, i), "w+") as out_file:
            for scaff_smi, dec_smi in sliced_mols_rdd.map(format_func).collect():
                out_file.write("{}\t{}\n".format(scaff_smi, dec_smi))


LOG = ul.get_logger("create_randomized_smiles")
if __name__ == "__main__":
    SPARK, SC = us.SparkSessionSingleton.get("create_randomized_smiles")
    main()
