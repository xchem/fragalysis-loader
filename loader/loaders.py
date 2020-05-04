import sys, json, os, glob
from django.contrib.auth.models import User
from viewer.models import Target, Protein, Molecule, Compound, Project
from hypothesis.models import (
    Vector3D,
    Vector,
    InteractionPoint,
    TargetResidue,
    ProteinResidue,
    Interaction,
)
from hypothesis.definitions import VectTypes, IntTypes
from hotspots.models import HotspotMap
from django.core.exceptions import ValidationError
from django.core.files import File
from rdkit import Chem
from rdkit.Chem import Descriptors
from scoring.models import MolGroup,MolAnnotation
from frag.alysis.run_clustering import run_lig_cluster
from loader.functions import sanitize_mol, get_path_or_none
from frag.network.decorate import get_3d_vects_for_mol
from loader.config import get_dict
import numpy as np
import pandas as pd
from rdkit.Geometry import Point3D


def add_target(title):
    """
    Add a target
    :param title: add a target by title
    :return: the created target
    """
    new_target = Target.objects.get_or_create(title=title)
    print('Target created = ' +str(Target.objects.get_or_create(title=title)[1]))
    return new_target[0]


def add_prot(pdb_file_path, code, target, mtz_path=None, map_path=None, bound_path=None):
    """
    Add a protein with a PDB, code and
    :param pdb_file_path: the PDB file path
    :param code: the unique code for this file
    :param target: the target to be linkede to
    :param mtz_path: the path to the MTZ file
    :param map_path: the path to the MAP file
    :return: the created protein
    """
    new_prot = Protein.objects.get_or_create(code=code, target_id=target)
    print('Protein created = ' + str(new_prot[1]))
    new_prot = new_prot[0]
    new_prot.apo_holo = True
    if pdb_file_path:
        new_prot.pdb_info.delete()
        new_prot.pdb_info.save(os.path.basename(pdb_file_path), File(open(pdb_file_path)))
    if mtz_path:
        new_prot.mtz_info.delete()
        new_prot.mtz_info.save(os.path.basename(mtz_path), File(open(mtz_path)))
    if map_path:
        new_prot.map_info.delete()
        new_prot.map_info.save(os.path.basename(map_path), File(open(map_path)))
    if bound_path:
        new_prot.bound_info.delete()
        new_prot.bound_info.save(os.path.basename(bound_path), File(open(bound_path)))
    new_prot.save()
    return new_prot


def add_projects_to_cmpd(new_comp, projects):
    """
    Add a project links to a compound
    :param new_comp: the Django compound to add them to
    :param projects:  the list Django projects to add
    :return: the compound with the added projects
    """
    [new_comp.project_id.add(x) for x in projects]
    new_comp.save()
    return new_comp


def calc_cpd(cpd_object, mol, projects):
    # Neutralise and desalt compound the compound
    sanitized_mol = sanitize_mol(mol)
    # Store the isomeric smiles
    smiles = Chem.MolToSmiles(sanitized_mol, isomericSmiles=True)
    # The inchi string is used for unique identification
    inchi = Chem.MolToInchi(sanitized_mol)
    # Now convert back to inchi to canonicalise
    tmp_mol = Chem.MolFromInchi(inchi)
    if tmp_mol is None:
        # If error in INNCHI READ -> NOT NECCESARILY A KILLER
        sys.stderr.write("INCHI ERROR: " + inchi)
    else:
        inchi = Chem.MolToInchi(tmp_mol)

    cpd_object.smiles = smiles
    if len(smiles) > Compound._meta.get_field("smiles").max_length:
        print("SMILES TOO LONG")
        return None
    if not len(inchi) > 255:
        cpd_object.inchi = inchi
    else:
        print("INCHI TOO LONG")
        return None
    m = sanitized_mol

    if m is None:
        sys.stderr.write("NONE MOLECULE PRODUCED\n" + smiles + "\n" + inchi)
        return None
    cpd_object.mol_log_p = Chem.Crippen.MolLogP(m)
    cpd_object.mol_wt = float(Chem.rdMolDescriptors.CalcExactMolWt(m))
    cpd_object.heavy_atom_count = Chem.Lipinski.HeavyAtomCount(m)
    cpd_object.heavy_atom_mol_wt = float(Descriptors.HeavyAtomMolWt(m))
    cpd_object.nhoh_count = Chem.Lipinski.NHOHCount(m)
    cpd_object.no_count = Chem.Lipinski.NOCount(m)
    cpd_object.num_h_acceptors = Chem.Lipinski.NumHAcceptors(m)
    cpd_object.num_h_donors = Chem.Lipinski.NumHDonors(m)
    cpd_object.num_het_atoms = Chem.Lipinski.NumHeteroatoms(m)
    cpd_object.num_rot_bonds = Chem.Lipinski.NumRotatableBonds(m)
    cpd_object.num_val_electrons = Descriptors.NumValenceElectrons(m)
    cpd_object.ring_count = Chem.Lipinski.RingCount(m)
    cpd_object.tpsa = Chem.rdMolDescriptors.CalcTPSA(m)
    # Validate that the compound is unique
    try:
        cpd_object.validate_unique()
        cpd_object.save()
        cpd_object = add_projects_to_cmpd(cpd_object, projects)
        return cpd_object
    except ValidationError:
        if not inchi:
            cpd_object.save()
        cpd_object = Compound.objects.get(inchi=inchi)
        cpd_object = add_projects_to_cmpd(cpd_object, projects)
        return cpd_object


def update_cpd(cpd_id, mol, projects):
    print(mol)
    cpd = cpd_id
    comp = calc_cpd(cpd, mol, projects)
    return comp


def add_comp(mol, projects, option=None, comp_id=None):
    """
    Function to add a new compound to the database given an RDKit molecule
    Takes an RDKit molecule.
    :param mol: the input RDKit molecule
    :param option: Option of LIG to return original smiles with the Compound object
    :param comp_id: the Django compound it relates to
    :return: a compound object for the RDKit molecule
    """

    # Now attribute all this meta-deta to the compound object
    new_comp = Compound()
    comp = calc_cpd(new_comp, mol, projects)
    return comp


def add_mol(mol_sd, prot, projects, lig_id="LIG", chaind_id="Z", occupancy=0.0):
    """
    Function to add a new Molecule to the database
    :param mol_sd: the SDMolBlock of the molecule
    :param prot: the protein it is associated to
    :param lig_id: the 3 letter ligand id
    :param chaind_id: the chain id
    :param occupancy: the occupancy
    :return: the created molecule
    """
    # create mol object from mol_sd
    rd_mol = Chem.MolFromMolFile(mol_sd)

    if rd_mol is None:
        return None

    # See if there is already a molecule with a compound
    old_mols = Molecule.objects.filter(prot_id=prot)

    print('OLD MOLS = ' + str(len(old_mols)))
    # If there's only one
    if len(old_mols)==1:
        # find the right id (if it exists)
        cpd_id = old_mols[0].cmpd_id
        if cpd_id:
            # update existing compound
            comp_ref = update_cpd(cpd_id, rd_mol, projects)
        else:
            # create new if there's no cpd already
            comp_ref = add_comp(rd_mol, projects)

    else:
        comp_ref = add_comp(rd_mol, projects)

    if comp_ref:
        new_mol = Molecule.objects.get_or_create(prot_id=prot, cmpd_id=comp_ref)[0]
        # Make a protein object by which it is related in the DB
        new_mol.sdf_info = Chem.MolToMolBlock(rd_mol)
        new_mol.smiles = Chem.MolToSmiles(rd_mol, isomericSmiles=True)
        # Find out how to add this information from Proasis
        new_mol.lig_id = lig_id
        new_mol.chain_id = chaind_id
        new_mol.occupancy = occupancy
        # Add this to the compound list -> make sure this passes in for the
        # correct molecule. I.e. if it fails where does it go???
        # Now link that compound back
        new_mol.cmpd_id = comp_ref
        new_mol.save()
        return new_mol
    else:
        return None


def parse_proasis(input_string):
    """
    Parse proasis contact strings
    :param input_string: the Proasis contact string to parse
    :return: a tuple of res_name, res_num, chain_id
    """
    return (
        input_string[:3].strip(),
        int(input_string[5:].strip()),
        input_string[3:5].strip(),
    )


def create_int(prot_res, mol, int_type, interaction):
    """
    Create a Django interaction object
    :param prot_res: the Django protein residue
    :param mol: the Django molecule
    :param int_type: the interaction type string
    :param interaction: the interaction dictionary
    :return: None
    """
    interation_point = InteractionPoint.objects.get_or_create(
        prot_res_id=prot_res,
        mol_id=mol,
        protein_atom_name=interaction["dstaname"],
        molecule_atom_name=interaction["srcaname"],
    )[0]
    Interaction.objects.get_or_create(
        interaction_version="PR",
        interaction_point=interation_point,
        interaction_type=int_type.get_int_conv("PR", interaction["contactType"]),
        distance=interaction["dis"],
        prot_smarts=interaction["dstType"],
        mol_smarts=interaction["srcType"],
    )


def add_contacts(input_data, target, prot, mol):
    """
    Add a series of Django contact objects
    :param input_data: the
    :param target: the data - either dict or list - of itneractions
    :param prot: the Django protein object
    :param mol: the Django molecule object
    :return: None
    """
    int_type = IntTypes()
    int_list = []
    if type(input_data) == dict:
        if "results" in input_data:
            int_list = input_data["results"]
    else:
        int_list = input_data
    for interaction in int_list:
        # Ignore Water mediasted hypothesis with Protein for now
        if interaction["hetmoltype"] == "WATER":
            continue
        res_name, res_num, chain_id = parse_proasis(interaction["dstrname"])
        targ_res = TargetResidue.objects.get_or_create(
            target_id=target, res_name=res_name, res_num=res_num, chain_id=chain_id
        )[0]
        prot_res = ProteinResidue.objects.get_or_create(
            targ_res_id=targ_res, prot_id=prot
        )[0]
        create_int(prot_res, mol, int_type, interaction)


def add_map(new_prot, new_target, map_path, map_type):
    """
    Add a Django map obect
    :param new_prot: the Django protein object
    :param new_target: the Django target object
    :param map_path: the path to the map file
    :param map_type: the two letter code signifyign the type of the map
    :return: the add Django map object
    """
    hotspot_map = HotspotMap.objects.get_or_create(
        map_type=map_type, target_id=new_target, prot_id=new_prot
    )[0]
    hotspot_map.map_info.save(os.path.basename(map_path), File(open(map_path)))
    return hotspot_map


def delete_users(project):
    """
    Refresh the users for a given project by deleting all of them.
    Redundant if using iSpyB.
    :param project: the project to remove users from.
    :return: None
    """
    for user_id in project.user_id.all():
        project.user_id.remove(user_id.pk)
    project.save()


def add_visits_or_proposal(target, file_path):
    """
    Add visits for a given target
    :param target: the target to add visits to
    :param file_path: the path to the file describing the available visits in space delimited format.
    :return: the Django projects created in this process
    """
    visits = [x.strip() for x in open(file_path).readlines() if x.strip()]
    projects = []
    for visit_line in visits:
        visit = visit_line.split()[0]
        project = Project.objects.get_or_create(title=visit)[0]
        projects.append(project)
        delete_users(project)
        target.project_id.add(project)
        for fedid in visit_line.split()[1:]:
            user = User.objects.get_or_create(username=fedid, password="")[0]
            project.user_id.add(user)
    target.save()
    return projects


def add_projects(new_target, dir_path, app):
    """
    Add proposals and visits as projects for a given target.
    :param new_target: the target being added
    :param dir_path: the path for where the PROPOSALS and VISITS files are held.
    :return: the projects added.
    """
    # Add the proposal information
    proposal_path = os.path.join(dir_path, "PROPOSALS")
    visit_path = os.path.join(dir_path, "VISITS")
    projects = []
    if os.path.isfile(proposal_path):
        projects.extend(add_visits_or_proposal(new_target, proposal_path))
    if os.path.isfile(visit_path):
        projects.extend(add_visits_or_proposal(new_target, visit_path))
    # remove_not_added(new_target, projects, app=app)
    return projects


def remove_not_added_projects(target, projects, app):
    """
    Remove any projects that have not been added this time around.
    Ensures the database updates, e.g. if projects or visits are added.
    :param target: the target added
    :param projects: the projects that have been added
    :return:
    """
    if app == 'fragspect':
        return None
    project_pks = [x.pk for x in projects]
    for project_id in target.project_id.all():
        if project_id.pk not in project_pks:
            target.project_id.remove(project_id.pk)
    target.save()


def remove_not_added(target, xtal_list, app):
    """
    Remove any crystals that have not been added this time around.
    Ensures the database updates, e.g. if someone nobody wants a given xtal.
    :param target: the target being considered
    :param xtal_list: a list of protein codes that have been added
    :return: None
    """
    if app == 'fragspect':
        return None
    all_prots = Protein.objects.filter(target_id=target)
    for prot in all_prots:
        if prot.code not in xtal_list:
            prot.delete()
    return None


def save_confidence(mol, file_path, annotation_type="ligand_confidence"):
    input_dict = json.load(open(file_path))
    val_store_dict = ["ligand_confidence_comment","refinement_outcome","ligand_confidence_int"]
    for val in val_store_dict:
        if val in input_dict:
            value = input_dict[val]
            if value:
                mol_annot = MolAnnotation.objects.get_or_create(mol_id=mol, annotation_type=annotation_type)[0]
                mol_annot.annotation_text = value
                mol_annot.save()
        else:
            print(val+ " not found in " + str(input_dict) + " for mol " + str(mol.prot_id.code))


def add_biomol_remark(search_path):
    biomol_remark = open(search_path + '/biomol.txt').readlines()
    print(biomol_remark)
    for f in glob.glob(search_path + '/*/*.pdb'):
        print(f)
        with open(f) as handle:
            switch = 0
            header_front, header_end = [], []
            pdb = []

            for line in handle:

                if line.startswith('ATOM'): switch = 1

                if line.startswith('HETATM'): switch = 2

                if switch == 0:
                    header_front.append(line)

                elif (switch == 2) and not line.startswith('HETATM'):
                    header_end.append(line)

                else:
                    pdb.append(line)
        full_file = ''.join(header_front) + ''.join(biomol_remark) + ''.join(pdb) + ''.join(header_end)
        with open(f, 'w') as w:
            w.write(full_file)


def load_from_dir(target_name, dir_path, app):
    """
    Load the data for a given target from a directory structure
    :param target_name: the string title of the target. This will uniquely identify it.
    :param dir_path: the path to the input data.
    :return: None
    """
    input_dict = get_dict()
    if os.path.isdir(dir_path):
        pass
    else:
        print("No data to add: " + target_name)
        return None
    new_target = add_target(target_name)
    projects = add_projects(new_target, dir_path, app=app)
    directories = sorted(os.listdir(dir_path))
    xtal_list = []
    for xtal in directories:
        if not os.path.isdir(os.path.join(dir_path,xtal)):
            continue
        print(xtal)
        xtal_list.append(xtal)
        new_path = os.path.join(dir_path, xtal)
        pdb_file_path = get_path_or_none(new_path, xtal, input_dict, "APO")
        bound_path = get_path_or_none(new_path, xtal, input_dict, "BOUND")
        mol_file_path = get_path_or_none(new_path, xtal, input_dict, "MOL")
        # using the pandda map for the target map file - for now
        map_path = get_path_or_none(new_path, xtal, input_dict, "PMAP")
        mtz_path = get_path_or_none(new_path, xtal, input_dict, "MTZ")
        # optional ones - contacts and hotspots
        contact_path = get_path_or_none(new_path, xtal, input_dict, "CONTACTS")
        ligand_confidence = get_path_or_none(new_path, xtal, input_dict, "CONFIDENCE")
        acc_path = get_path_or_none(new_path, xtal, input_dict, "ACC")
        don_path = get_path_or_none(new_path, xtal, input_dict, "DON")
        lip_path = get_path_or_none(new_path, xtal, input_dict, "LIP")
        if pdb_file_path and mol_file_path:
            if os.path.isfile(pdb_file_path) and os.path.isfile(mol_file_path):
                new_prot = add_prot(
                    pdb_file_path, xtal, new_target, mtz_path=mtz_path, map_path=map_path, bound_path=bound_path
                )
                new_mol = add_mol(mol_file_path, new_prot, projects)
                if not new_mol:
                    print("NONE MOL: " + xtal)
                else:
                    if contact_path:
                        try:
                            add_contacts(
                                json.load(open(contact_path)), new_target, new_prot, new_mol
                            )
                        except ValueError:
                            print("Error parsing: " + contact_path)
                    else:
                        print("Skipping contacts - " + xtal)
                    if ligand_confidence:
                        save_confidence(new_mol,ligand_confidence)
                    else:
                        print("Skipping confidence - " + xtal)
                    if acc_path:
                        add_map(new_prot, new_target, acc_path, "AC")
                    if don_path:
                        add_map(new_prot, new_target, don_path, "DO")
                    if lip_path:
                        add_map(new_prot, new_target, lip_path, "AP")
        elif bound_path and map_path:
            if os.path.isfile(bound_path) and os.path.isfile(map_path):
                new_prot = add_prot(
                    pdb_file_path, xtal, new_target, mtz_path=mtz_path, map_path=map_path, bound_path=bound_path
                )
        else:
            print("File not found: " + xtal)
    remove_not_added(new_target, xtal_list, app=app)


def create_vect_3d(mol, new_vect, vect_ind, vector):
    """
    Generate the 3D synthesis vectors for a given molecule
    :param mol: the Django molecule object
    :param new_vect: the Django 2d vector object
    :param vect_ind: the index of the vector - since the same 2D vector
    can be different in 3D
    :param vector: the vector coordinates - a 2*3 list of lists.
    :return: None
    """
    if vector:
        new_vect3d = Vector3D.objects.get_or_create(
            mol_id=mol, vector_id=new_vect, number=vect_ind
        )[0]
        # The start position
        new_vect3d.start_x = float(vector[0][0])
        new_vect3d.start_y = float(vector[0][1])
        new_vect3d.start_z = float(vector[0][2])
        # The end position
        new_vect3d.end_x = float(vector[1][0])
        new_vect3d.end_y = float(vector[1][1])
        new_vect3d.end_z = float(vector[1][2])
        new_vect3d.save()


def get_vectors(mols):
    """
    Get the vectors for a given molecule
    :param mols: the Django molecules to get them from
    :param target: the Django target to record them from
    :return: None
    """
    vect_types = VectTypes()
    for mol in mols:
        if "." in mol.smiles:
            print("SKIPPING - FRAGMENT: " + str(mol.pk)) + " " + str(mol.smiles)
            continue
        vectors = get_3d_vects_for_mol(mol.sdf_info)
        for vect_type in vectors:
            vect_choice = vect_types.translate_vect_types(vect_type)
            for vector in vectors[vect_type]:
                spl_vect = vector.split("__")
                smiles = spl_vect[0]
                if len(spl_vect) > 1:
                    vect_ind = int(spl_vect[1])
                else:
                    vect_ind = 0
                new_vect = Vector.objects.get_or_create(
                    smiles=smiles, cmpd_id=mol.cmpd_id, type=vect_choice
                )[0]
                create_vect_3d(mol, new_vect, vect_ind, vectors[vect_type][vector])


def cluster_mols(rd_mols, mols, target):
    """
    Cluster a series of 3D molecules
    :param rd_mols: the RDKit molecules to cluster
    :param mols:  the Django moleculs they refer to
    :param target:  the Django target it refers to
    :return: None
    """
    id_mols = [x.pk for x in mols]
    out_data = run_lig_cluster(rd_mols, id_mols)
    for clust_type in out_data:
        for cluster in out_data[clust_type]:
            # look for molgroup with same coords - need to implement tolerance?
            mol_group = search_for_molgroup_by_coords(coords=[out_data[clust_type][cluster]["centre_of_mass"][0],
                                                              out_data[clust_type][cluster]["centre_of_mass"][1],
                                                              out_data[clust_type][cluster]["centre_of_mass"][2]],
                                                      target=target.title)
            if not mol_group:
                mol_group = MolGroup()
            if clust_type != "c_of_m":
                mol_group.group_type = "PC"
            else:
                mol_group.group_type = "MC"
            mol_group.target_id = target
            mol_group.x_com = out_data[clust_type][cluster]["centre_of_mass"][0]
            mol_group.y_com = out_data[clust_type][cluster]["centre_of_mass"][1]
            mol_group.z_com = out_data[clust_type][cluster]["centre_of_mass"][2]
            mol_group.description = clust_type
            mol_group.save()
            for mol_id in out_data[clust_type][cluster]["mol_ids"]:
                if mol_id not in [a['id'] for a in mol_group.mol_id.values()]:
                    this_mol = Molecule.objects.get(id=mol_id)
                    mol_group.mol_id.add(this_mol)


def centre_of_points(list_of_points):
    cp = np.average(list_of_points, axis=0)
    return cp


def centre_of_mass(mol):
    numatoms = mol.GetNumAtoms()
    conf = mol.GetConformer()
    if not conf.Is3D():
        return 0
    # get coordinate of each atoms
    pts = np.array([list(conf.GetAtomPosition(atmidx)) for atmidx in range(numatoms)])
    atoms = [atom for atom in mol.GetAtoms()]
    mass = Descriptors.MolWt(mol)
    # get center of mass
    center_of_mass = np.array(np.sum(atoms[i].GetMass() * pts[i] for i in range(numatoms))) / mass
    return center_of_mass


def process_site(rd_mols):
    coms = [centre_of_mass(mol) for mol in rd_mols]
    centre = centre_of_points(coms)
    print('CENTRE: ' + str(centre))
    return centre


def get_coord_limits(coord):
    lower_limit = float('.'.join([str(coord).split('.')[0], str(coord).split('.')[1][:2]]))
    if lower_limit > 0:
        upper_limit = lower_limit + 0.01
    else:
        tmp = lower_limit - 0.01
        upper_limit = lower_limit
        lower_limit = tmp
    return lower_limit, upper_limit


def search_for_molgroup_by_coords(coords, target):
    x = coords[0]
    y = coords[1]
    z = coords[2]

    limit_list = []

    for coord in x, y, z:
        lower, upper = get_coord_limits(coord)
        limit_list.append([lower, upper])

    search = MolGroup.objects.filter(target_id__title=target, x_com__gte=limit_list[0][0], x_com__lte=limit_list[0][1],
                                     y_com__gte=limit_list[1][0], y_com__lte=limit_list[1][1],
                                     z_com__gte=limit_list[2][0],
                                     z_com__lte=limit_list[2][1])

    if len(search) == 1:
        mol_group = search[0]
    else:
        return None

    return mol_group


def search_for_molgroup_by_description(description, target):
    search = MolGroup.objects.filter(target_id__title=target, description=description)
    print(str('matching_sites = ')+str(len(search)))
    if len(search)==1:
        mol_group = search[0]

    elif len(search) >1:
        for molgroup in search:
            molgroup.delete()
        return None
    else:
        return None

    return mol_group


def analyse_mols(mols, target, specified_site=False, site_description=None):
    """
    Analyse a list of molecules for a given target
    :param mols: the Django molecules to analyse
    :param target: the Django target
    :return: None
    """
    rd_mols = [Chem.MolFromMolBlock(x.sdf_info) for x in mols]
    if not specified_site:

        cluster_mols(rd_mols, mols, target)

    else:

        centre = process_site(rd_mols)

        # look for molgroup with same target and description
        mol_group = search_for_molgroup_by_description(target=target.title, description=site_description)

        if not mol_group:
            mol_group = MolGroup()
        mol_group.group_type = "MC"
        mol_group.target_id = target
        mol_group.x_com = centre[0]
        mol_group.y_com = centre[1]
        mol_group.z_com = centre[2]
        mol_group.description = site_description
        mol_group.save()

        ids = [m.id for m in mols]

        print([a['id'] for a in mol_group.mol_id.values()])

        for mol_id in ids:
            if mol_id not in [a['id'] for a in mol_group.mol_id.values()]:
                print(mol_id)
                this_mol = Molecule.objects.get(id=mol_id)
                mol_group.mol_id.add(this_mol)

    get_vectors(mols)



def analyse_target(target_name, target_path):
    """
    Analyse all the molecules for a particular target
    :param target_name: the name of the target
    :return: None
    """
    target = Target.objects.get(title=target_name)
    mols = list(Molecule.objects.filter(prot_id__target_id=target))
    print("Analysing " + str(len(mols)) + " molecules for " + target_name)
    # Delete the old ones for this target - don't delete! UPDATE...
    # MolGroup.objects.filter(group_type="PC", target_id=target).delete()
    # MolGroup.objects.filter(group_type="MC", target_id=target).delete()
    if os.path.isfile(os.path.join(target_path, 'hits_ids.csv')) and os.path.isfile(
            os.path.join(target_path, 'sites.csv')):

        hits_sites = pd.read_csv(os.path.join(target_path, 'hits_ids.csv'))
        sites = pd.read_csv(os.path.join(target_path, 'sites.csv'))

        for i, row in sites.iterrows():
            description = row['site']
            print('Processing user input site: ' + description)
            hit_ids = list(hits_sites['crystal_id'][hits_sites['site_number'] == i])
            print('HIT IDS: ' + str(hit_ids))
            mols = list(Molecule.objects.filter(prot_id__target_id=target, prot_id__code__in=hit_ids))
            analyse_mols(mols=mols, target=target, specified_site=True, site_description=description)

    else:
        analyse_mols(mols=mols, target=target)

def write_lig_pdb(pdb):
    f = open('lig.pdb', 'w')
    lig_atoms = []
    link_info = None
    for line in pdb:
        if line.startswith('LINK'):
            f.write(line)
            link_info = [line[13]]+[line[42:46].strip()]+line[47:62].strip().split(' ')
        if line[17:20] == 'LIG' and line.startswith('HETATM'):
            lig_atoms.append(line[7:11].strip())
            f.write(line)
        if line.startswith('CONECT'):
            yes = False
            for num in line.split(' '):
                if num in lig_atoms:
                    yes = True
            if yes == True:
                f.write(line)
    f.close()
    return link_info


def get_3d_distance(coord_a, coord_b):
    sum_ = (sum([(float(coord_a[i])-float(coord_b[i]))**2 for i in range(3)]))
    return np.sqrt(sum_)


def process_covalent(directory):
    link_atoms = {'SG': 16, 'O': 8, 'N': 7, 'C':12, 'C5':12}
    print([x[0] for x in os.walk(directory)])
    for f in [x[0] for x in os.walk(directory)]:
        print(str(f) + '/*_bound.pdb')
        print(glob.glob(str(f) + '/*_bound.pdb'))
        if glob.glob(str(f) + '/*_bound.pdb'):
            bound_pdb = glob.glob(str(f) + '/*_bound.pdb')[0]
            mol_file = glob.glob(str(f) + '/*.mol')[0]
            pdb = open(bound_pdb, 'r').readlines()

            link_info = write_lig_pdb(pdb)
            if link_info is not None:
                m2 = Chem.MolFromPDBFile('lig.pdb')

                for line in pdb:
                    if line[13:17].strip() == link_info[1] and line[17:20] == link_info[2] and line[20:23].strip() == \
                            link_info[3] and line[23:27].strip() == link_info[4]:
                        # res_coords = [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()]
                        try:
                            res_coords = [float(line[31:39].strip()), float(line[39:47].strip()),
                                          float(line[47:55].strip())]
                        except ValueError:
                            pass

                edmol = Chem.EditableMol(m2)
                new_mol = edmol.AddAtom(Chem.Atom(link_atoms[link_info[1]]))
                Chem.MolToPDBFile(edmol.GetMol(), 'edlig.pdb')
                edpdb = open('edlig.pdb', 'r').readlines()

                distances = {}
                for line in edpdb:
                    if line.startswith('HETATM'):
                        distances[get_3d_distance(res_coords,
                                                  [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()])] = line[
                                                                                                                      7:11].strip()
                    if [line[31:39].strip(), line[39:47].strip(), line[47:55].strip()] == ['0.000', '0.000', '0.000']:
                        new_idx = line[7:11].strip()

                edmol.AddBond(int(distances[min(distances)]) - 1, int(new_idx) - 1, Chem.BondType.SINGLE)

                new_mol = edmol.GetMol()
                conf = new_mol.GetConformer()
                res_coords = tuple([float(i) for i in res_coords])
                conf.SetAtomPosition(int(new_idx) - 1, Point3D(res_coords[0], res_coords[1], res_coords[2]))
                Chem.MolToMolFile(new_mol, mol_file)


def process_target(prefix, target_name, app):
    """
    Process the full target.
    :param prefix:
    :param target_name:
    :return:
    """
    target_path = os.path.join(prefix, target_name)
    print('TARGET_PATH: ' + target_path)
    process_covalent(target_path)
    if os.path.isfile(os.path.join(target_path, 'biomol.txt')):
        print('ADDING BIOMOL REMARK')
        add_biomol_remark(search_path=target_path)
    load_from_dir(target_name, target_path, app=app)
    # Check for new data
    new_data_file = os.path.join(target_path, "NEW_DATA")
    if os.path.isfile(new_data_file) and app == 'fragspect':
        os.remove(new_data_file)

    if os.path.isfile(new_data_file):
        print("Analysing target: " + target_name)
        analyse_target(target_name, target_path)
        os.remove(new_data_file)
    else:
        print("NEW_DATA not found for " + target_path)
