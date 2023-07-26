#! /bin/bash 

set -e

cargo build

rm book/src/advanced_usage.md
cp advanced_usage.template book/src/advanced_usage.md

doc_file=book/src/advanced_usage.md

echo "" >> ${doc_file}
echo "" >> ${doc_file}
echo "\`\`\`text" >> ${doc_file}
./target/debug/modkit --help >> ${doc_file}
echo "\`\`\`" >> ${doc_file}

subcommands=("pileup" "adjust-mods" "update-tags" "sample-probs" "summary" "motif-bed" "call-mods" "extract" "repair")
for cmd in "${subcommands[@]}"
do
  echo "" >> ${doc_file}
  echo "## ${cmd}" >> ${doc_file}
  echo "\`\`\`text" >> ${doc_file}
  ./target/debug/modkit $cmd --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

