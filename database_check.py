from viewer.models import Molecule, Compound, Protein, Target, Project
import sys,os

class RunDbCheck(object):

    def run_check(self):
        """

        :return:
        """
        self.project = Project.objects.create(title="DUMMY_PROJECT")
        self.target = Target.objects.create(title="DUMMY_TARGET")
        self.target.project_id.add(self.project)
        self.target.save()
        self.cmpd = Compound.objects.create(
            inchi="DUM_INCH",
            smiles=self.mol_smi,
            mol_log_p=0.1,
            mol_wt=0.2,
            tpsa=0.3,
            heavy_atom_count=1,
            heavy_atom_mol_wt=2,
            nhoh_count=3,
            no_count=4,
            num_h_acceptors=5,
            num_h_donors=6,
            num_het_atoms=7,
            num_rot_bonds=8,
            num_val_electrons=9,
            ring_count=10,
        )
        self.cmpd_two = Compound.objects.create(
            inchi="DUM_INCH_TWO",
            smiles=self.mol_smi_two,
            mol_log_p=0.1,
            mol_wt=0.2,
            tpsa=0.3,
            heavy_atom_count=1,
            heavy_atom_mol_wt=2,
            nhoh_count=3,
            no_count=4,
            num_h_acceptors=5,
            num_h_donors=6,
            num_het_atoms=7,
            num_rot_bonds=8,
            num_val_electrons=9,
            ring_count=10,
        )
        self.protein = Protein.objects.create(
            code="DUMM", target_id=self.target, pdb_info="my_pdb.pdb"
        )
        self.protein_two = Protein.objects.create(
            code="DUMM_TWO", target_id=self.target, pdb_info="my_pdb.pdb"
        )
        self.mol = Molecule.objects.create(
            smiles="O=C(c1ccc2c(c1)OCO2)N1CCCCCC1",
            lig_id="DUM",
            chain_id="C",
            sdf_info=self.mol_sd_str,
            rscc=0.1,
            occupancy=0.2,
            x_com=0.3,
            y_com=0.4,
            z_com=0.5,
            rmsd=0.6,
            prot_id=self.protein,
            cmpd_id=self.cmpd,
        )
        self.dj_mol_two = Molecule.objects.create(
            smiles=self.mol_smi_two,
            lig_id="DUM",
            chain_id="C",
            sdf_info=self.mol_sd_two_str,
            rscc=0.1,
            occupancy=0.2,
            x_com=0.3,
            y_com=0.4,
            z_com=0.5,
            rmsd=0.6,
            prot_id=self.protein_two,
            cmpd_id=self.cmpd_two,
        )

    def run_test(self):
        """
        Test that the data has been added appropriately.
        :return:
        """
        Project.objects.get(title="DUMMY_PROJECT")
        Target.objects.get(title="DUMMY_TARGET")

    def run_delete(self):
        """
        Delete the added data
        :return:
        """
        # Delete the data just added
        self.project.delete()
        self.target.delete()

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django

    django.setup()

    run_db_check = RunDbCheck()
    run_db_check.run_check()
    run_db_check.run_test()
    run_db_check.run_delete()