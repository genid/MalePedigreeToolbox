"""Testing of commands with inputs as a sanity check for command functionality"""

from unittest import TestCase

from MalePedigreeToolbox.tests.testing_utility import *


class Test(TestCase):

    def setUp(self) -> None:
        create_temp_out()

    def tearDown(self) -> None:
        clean_temp_out()
        # pass

    @classmethod
    def tearDownClass(cls):
        clean_log_files()
        # pass

    def test_distance_command(self):
        output_dir = TEMP_OUT_DIR / 'test_outdir'
        command = f'mpt -f -ll silent distances -t "{TEST_FILE_DIR / "small_TGFs"}" -o "{output_dir}"'
        run_command(command)
        self.assertTrue(confirm_files_exist(output_dir / "distances.csv"))
        self.assertTrue(confirm_lines_equal(output_dir / "distances.csv", TEST_FILE_DIR / 'expected_distance_out.csv'))

    def test_mutation_diff_command(self):
        output_dir = TEMP_OUT_DIR / 'test_outdir'
        command = f'mpt -f -ll silent pairwise_mutation -af "{TEST_FILE_DIR / "alleles_small.csv"}" -df ' \
                  f'"{TEST_FILE_DIR / "distances_mut_diff_test.csv"}" -o "{output_dir}"'
        run_command(command)
        self.assertTrue(confirm_files_exist(output_dir / "full_out.csv", output_dir / "summary_out.csv",
                                            output_dir / "differentiation_out.csv"))
        self.assertTrue(confirm_columns_equal(output_dir / "full_out.csv", TEST_FILE_DIR / "expected_fo.csv", [0]))
        self.assertTrue(confirm_columns_equal(output_dir / "summary_out.csv", TEST_FILE_DIR / "expected_so.csv", [0]))
        self.assertTrue(confirm_lines_equal(output_dir / "differentiation_out.csv", TEST_FILE_DIR / "expected_do.csv"))

    def test_infer_pedigree_mutations(self):
        output_dir = TEMP_OUT_DIR / 'test_outdir'
        command = f'mpt -f -ll silent pedigree_mutation -t "{TEST_FILE_DIR / "infer_pedigree_tgfs"}" -af ' \
                  f'"{TEST_FILE_DIR / "infer_ped_mut_alleles.csv"}"' \
                  f' -o "{output_dir}"'
        run_command(command)
        for name, expected in [["chained_mutations", ["6", "2"]],
                               ["correct_split", ["10", "7"]],
                               ["correct_start_allele", ["4", "4"]],
                               ["one_mutation", ["5", "3"]],
                               ["simple_no_mutation", ["1", "0"]]]:
            out_folder = output_dir / name
            if name == "simple_no_mutation":
                self.assertTrue(confirm_files_exist(out_folder / f"{name}_mutations.csv",
                                                    out_folder / f"pedigree_{name}_all_marker_edge_info.csv",
                                                    out_folder / f"pedigree_{name}_all_markers.pdf"))
            else:
                self.assertTrue(confirm_files_exist(out_folder / f"{name}_mutations.csv",
                                                    out_folder / f"pedigree_{name}_all_marker_edge_info.csv",
                                                    out_folder / f"pedigree_{name}_all_markers.pdf",
                                                    out_folder / f"pedigree_{name}_marker_marker1.pdf"))
            with open(out_folder / f"{name}_mutations.csv") as f:
                lines = f.readlines()
                values = lines[1].split(",")
                self.assertEqual(values[1], expected[0])
                self.assertEqual(values[2], expected[1])

    def test_predict_pedigrees(self):
        output_dir = TEMP_OUT_DIR / 'test_outdir'
        command = f'mpt -f -ll silent dendrograms -c opt -fm "{TEST_FILE_DIR / "expected_fo.csv"}" ' \
                  f'-mr "{TEST_FILE_DIR / "marker_rates.csv"}" -o "{output_dir}"'
        run_command(command)
        self.assertTrue(confirm_files_exist(output_dir / "Draulans" / "Draulans_dendogram_clusters.txt"),
                        "no clusters file")
        self.assertTrue(confirm_files_exist(output_dir / "Draulans" / "Draulans_predicted_dendogram.png"),
                        "no dendrogram image file")
        self.assertTrue(confirm_lines_equal(output_dir / "Draulans" / "Draulans_dendogram_clusters.txt",
                                            TEST_FILE_DIR / "expected_dendrogram_clusters.txt"))

    def test_all(self):
        output_dir = TEMP_OUT_DIR / 'test_outdir'
        command = f'mpt -f -ll silent all -c opt -t "{TEST_FILE_DIR / "infer_pedigree_tgfs"}" -af ' \
                  f'"{TEST_FILE_DIR / "infer_ped_mut_alleles.csv"}" -o "{output_dir}"'
        run_command(command)
        self.assertTrue(confirm_files_exist(output_dir / "differentiation_out.csv", output_dir / "distances.csv",
                                            output_dir / "full_out.csv", output_dir / "summary_out.csv",
                                            output_dir / "total_mutations.csv"))
        for name, expected in [["chained_mutations", ["6", "2"]],
                               ["correct_split", ["10", "7"]],
                               ["correct_start_allele", ["4", "4"]],
                               ["one_mutation", ["5", "3"]],
                               ["simple_no_mutation", ["1", "0"]]]:
            out_folder = output_dir / name
            if name == "simple_no_mutation":
                self.assertTrue(confirm_files_exist(out_folder / f"{name}_mutations.csv",
                                                    out_folder / f"pedigree_{name}_all_marker_edge_info.csv",
                                                    out_folder / f"pedigree_{name}_all_markers.pdf",
                                                    out_folder / f"{name}_dendogram_clusters.txt",
                                                    out_folder / f"{name}_predicted_dendogram.png"))
            else:
                self.assertTrue(confirm_files_exist(out_folder / f"{name}_mutations.csv",
                                                    out_folder / f"pedigree_{name}_all_marker_edge_info.csv",
                                                    out_folder / f"pedigree_{name}_all_markers.pdf",
                                                    out_folder / f"pedigree_{name}_marker_marker1.pdf",
                                                    out_folder / f"{name}_dendogram_clusters.txt",
                                                    out_folder / f"{name}_predicted_dendogram.png"))
