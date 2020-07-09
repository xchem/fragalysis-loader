# Fragalysis Mounted Volume Data Loader Image
An extension of the `xchem/fragalysis-loader` image that first
copies files from a volume mounted at `/mounted-media` before calling
the loader's `run_loader.sh` script.
