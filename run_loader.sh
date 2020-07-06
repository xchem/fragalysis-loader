#!/bin/bash
export DJANGO_SETTINGS_MODULE=fragalysis.settings
echo "Running loader..."
python loader.py
script="
from viewer.models import Target,Project,Protein
from django.core.files import File
example_project = Project.objects.get_or_create(title='private_dummy_project')[0]
example_target = Target.objects.get_or_create(title='private_dummy_target')[0]
example_target.project_id.add(example_project)
example_target.save()
pdb_file_path = '/code/dummy_pdb.pdb'
out_f = open(pdb_file_path,'w')
out_f.write('IF YOURE READING THIS - STOP')
out_f.close()
new_prot = Protein.objects.get_or_create(code='dummy_protein_private', target_id=example_target)[0]
new_prot.apo_holo = True
new_prot.pdb_info.save(os.path.basename(pdb_file_path), File(open(pdb_file_path)))
print('Generated private data')
"
printf "$script" | python manage.py shell