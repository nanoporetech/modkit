#! /bin/bash 

set -e

cargo build

rm advanced_usage.md
cp advanced_usage.template advanced_usage.md

doc_file=advanced_usage.md

echo "" >> ${doc_file}
echo "" >> ${doc_file}
echo "\`\`\`bash" >> ${doc_file}
./target/debug/modkit --help >> ${doc_file}
echo "\`\`\`" >> ${doc_file}

subcommands=("pileup" "adjust-mods" "update-tags" "sample-probs" "summary" "motif-bed")
for cmd in "${subcommands[@]}"
do
  echo "" >> ${doc_file}
  echo "## ${cmd}" >> ${doc_file}
  echo "\`\`\`bash" >> ${doc_file}
  ./target/debug/modkit $cmd --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

