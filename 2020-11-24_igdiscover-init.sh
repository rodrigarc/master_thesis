#!/bin/bash -l

foo() {
echo "$f"
igdiscover init --db "/Users/karlor/Box\ Sync/RSV\ NGS/Database\ comparison/databases_uppmax/$DB" --single-reads "/Users/karlor/Box\ Sync/Rodrigo\ \&\ Klara\ \&\ Sebastian/Quality\ reports/Heavy_chains/fasta/combined/$f" "/Users/karlor/Box\ Sync/RSV\ NGS/RSV/germline_dbs/results/2020-11-23/IMGTrm/${f::3}"}
DB=$1

for f in E*; do foo $f; done
