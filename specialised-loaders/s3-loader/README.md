# Fragalysis S3 Data Loader Image
An extension of the `xchem/fragalysis-loader` image that first
copies files form an S3 bucket before calling the loader's
`run_loader.sh` script.

Due to the load times for large numbers of small files the media data
that resides on S3 is stored compressed in a single file (`media.tar.gz`)
and downloaded and decompressed, in-situ, by the loader.

If you have data in the directory `2020-09-15T16` compress it like this: -

    $ cd 2020-09-15T16
    $ tar -zcvf ../media.tar.gz .
    
...and then upload the `media.tar.gz` to the bucket path `2020-09-15T16`.

---
