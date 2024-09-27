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

subcommands=(
  "pileup"
  "adjust-mods"
  "update-tags"
  "sample-probs"
  "summary"
  "call-mods"
  "repair"
  "validate"
  "pileup-hemi"
  "entropy"
  "localize"
  "stats"
)

for cmd in "${subcommands[@]}"
do
  echo "" >> ${doc_file}
  echo "## ${cmd}" >> ${doc_file}
  echo "\`\`\`text" >> ${doc_file}
  ./target/debug/modkit $cmd --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

for subcommand in "full" "calls"; do
  echo "" >> ${doc_file}
  echo "## extract ${subcommand}" >> ${doc_file}
  echo "\`\`\`text" >> ${doc_file}
  ./target/debug/modkit extract $subcommand --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

for subcommand in "bed" "search" "evaluate" "refine"; do
  echo "" >> ${doc_file}
  echo "## motif ${subcommand}" >> ${doc_file}
  echo "\`\`\`text" >> ${doc_file}
  ./target/debug/modkit motif $subcommand --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

for subcommand in "pair" "multi"
do
  echo "" >> ${doc_file}
  echo "## dmr ${subcommand}" >> ${doc_file}
  echo "\`\`\`text" >> ${doc_file}
  ./target/debug/modkit dmr $subcommand --help >> ${doc_file}
  echo "\`\`\`" >> ${doc_file}
done

