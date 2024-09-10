#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from math import ceil
import tempfile
import os, shutil
from os.path import join, abspath
import subprocess
import itertools
from scipy.sparse import csr_array
import pandas as pd


@dataclass
class GfaseMaxcutSolver:
    solve_executable: str | None = "solve_maxcut"
    solver_image: str | None = None

    def _write_contacts(self, allele_matrix: csr_array, output_path: str) -> None:
        if not output_path.endswith(".csv"):
            raise ValueError()
        contacts = []
        assert allele_matrix.shape is not None
        for i in range(allele_matrix.shape[0]):
            nonzero_columns = allele_matrix[[i], :].nonzero()[1]  # type: ignore
            weight = ceil(1000 / len(nonzero_columns) ** 2)
            for j1, j2 in itertools.combinations(nonzero_columns, 2):
                a1, a2 = allele_matrix[i, j1] - 1, allele_matrix[i, j2] - 1  # type: ignore
                contacts.append((j1, a1, j2, a2, weight))

        PREFIX = "PR.0."
        with open(output_path, "wt") as f:
            f.write("name_a,name_b,weight\n")
            for j1, a1, j2, a2, weight in contacts:
                line = f"{PREFIX}{j1}.{a1},{PREFIX}{j2}.{a2},{weight}\n"
                f.write(line)

    def _write_ids(self, variant_count: int, output_path: str) -> None:
        if not output_path.endswith(".txt"):
            raise ValueError()
        PREFIX = "PR.0."
        with open(output_path, "wt") as f:
            for j in range(variant_count):
                f.write(f"{j * 2},{PREFIX}{j}.0\n")
                f.write(f"{j * 2 + 1},{PREFIX}{j}.1\n")

    def _load_haplotypes(self, input_path: str):
        if not input_path.endswith(".csv"):
            raise ValueError()
        df = pd.read_csv(input_path)
        haplotype_dict = {}
        for side, nodes in zip(df["side"], df["nodes"]):
            for node in nodes.split(" "):
                if not node:
                    continue
                fields = node.split(".")
                j = int(fields[2])
                a = int(fields[3])
                haplotype_dict[(j, a)] = side
        return haplotype_dict

    def solve(
        self,
        allele_matrix: csr_array,
        *,
        core_iterations: int = 200,
        sample_size: int = 50,
        n_rounds: int = 3,
        threads: int = 8,
        temp_dir: str | None = None,
        keep_intermediates: bool = False,
        verbose: bool = False,
    ):
        assert allele_matrix.shape is not None
        if temp_dir is not None:
            os.makedirs(temp_dir, exist_ok=True)
        workdir = tempfile.TemporaryDirectory(
            prefix="GfaseMaxcutSolver_", delete=False, dir=temp_dir
        ).name

        contacts_path = join(workdir, "contacts.csv")
        ids_path = join(workdir, "ids.txt")
        output_dir = join(workdir, "output")
        os.makedirs(output_dir)
        if self.solver_image:
            command = [
                "docker",
                "run",
                "-v",
                f"{contacts_path}:/data/contacts.csv",
                "-v",
                f"{ids_path}:/data/ids.txt",
                "-v",
                f"{output_dir}:/data/output",
                self.solver_image,
                "solve_maxcut",
                "-i",
                "/data/ids.txt",
                "-g",
                "/data/contacts.csv",
                "-o",
                "/data/output",
                "-c",
                core_iterations,
                "-s",
                sample_size,
                "-r",
                n_rounds,
                "-t",
                threads,
            ]
        elif self.solve_executable:
            command = [
                self.solve_executable,
                "-i",
                ids_path,
                "-g",
                contacts_path,
                "-o",
                output_dir,
                "-c",
                core_iterations,
                "-s",
                sample_size,
                "-r",
                n_rounds,
                "-t",
                threads,
            ]
        else:
            raise ValueError()
        command = [str(x) for x in command]

        try:
            self._write_ids(allele_matrix.shape[1], ids_path)
            self._write_contacts(allele_matrix, contacts_path)
            p = subprocess.run(
                command,
                capture_output=verbose,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            output_file = join(output_dir, "components_final.csv")
            if not os.path.isfile(output_file):
                print(p.stdout)
                print(p.stderr)
                raise RuntimeError("Failed to run GFAse")
            haplotype_dict = self._load_haplotypes(output_file)
        finally:
            if not keep_intermediates:
                shutil.rmtree(workdir)
        return haplotype_dict
