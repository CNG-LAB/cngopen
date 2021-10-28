#!/bin/bash

List=$(seq 1 ${2})
for i in $List
do
  echo node_${i} >> ${3}${1}.header
done

