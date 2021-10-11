#!/usr/bin/env bash

while read -r mdir op mfile
do
    if [ ! -f $1 ]
    then
        break
    fi

    touch $1
done < <(exec -a inotify_touch inotifywait -m -r -e modify -e delete -e delete_self -e move -e move_self . )
