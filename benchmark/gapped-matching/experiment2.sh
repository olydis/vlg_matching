#!/bin/bash

col=sources
prefix_size='1048576 20971520 209715200'
tmp_col='tmp.raw'

for size in ${prefix_size}; do
    echo "collections/$col/$col.raw"
    head -c ${size} collections/$col/$col.raw > ${tmp_col}
    ls -la ${tmp_col}
    echo "collections/$col-$size/"
    build/create_collection.x -i ${tmp_col} -c collections/$col-$size
    rm ${tmp_col}
done


# TODO: execute experiment
