import os
import glob


def get_target_list(base_path, ln):
    target_path = os.path.join(base_path, ln)
    if os.path.isfile(target_path):
        return [x.strip() for x in open(target_path).read().split(" ") if x.strip()]
    elif ln in os.environ:
        return os.environ[ln].split(",")
    else:
        return ["MURD", "HAO1A", "smTGR", "PTP1B"]


if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "fragalysis.settings")
    import django

    django.setup()
    from loader.loaders import process_target
    from loader.compound_set_loaders import process_compound_set

    prefix = "/code/media/NEW_DATA/"

    fragalysis_list = 'TARGET_LIST'
    fragspect_list = 'FRAGSPECT_LIST'

    if os.path.isfile(os.path.join(prefix, fragalysis_list)):
        list_name = fragalysis_list
        app = 'fragalysis'

    elif os.path.isfile(os.path.join(prefix, fragspect_list)):
        list_name = fragspect_list
        app = 'fragspect'

    else:
        raise Exception('No target list for Fragalysis or Fragspect found!')

    targets_to_load = get_target_list(prefix, ln=list_name)

    for target_name in targets_to_load:
        process_target(prefix, target_name, app=app)

    for target_name in targets_to_load:
        print(os.path.join(prefix, target_name, 'compound_sets/compound-set*.sdf'))
        compound_sets = glob.glob(os.path.join(prefix, target_name, 'compound_sets/compound-set*.sdf'))
        for cset in compound_sets:
            process_compound_set(target=target_name, filename=cset)
