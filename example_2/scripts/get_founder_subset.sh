#!/bin/bash

touch $1.subsampled

while read line; do
	grep "$line" $1 >> $1.subsampled
done < $2
