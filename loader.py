import os


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

    prefix = "/code/media/NEW_DATA/"

    fragalysis_list = 'TARGET_LIST'
    fragspect_list = 'FRAGSPECT_LIST'

    if os.path.isfile(os.path.join(prefix, fragalysis_list)):
        list_name = fragalysis_list
    elif os.path.isfile(os.path.join(prefix, fragspect_list)):
        list_name = fragspect_list

    targets_to_load = get_target_list(prefix, ln=list_name)

    for target_name in targets_to_load:
        process_target(prefix, target_name)

