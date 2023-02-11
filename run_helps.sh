#! /bin/bash 

set -e

cargo build

rm README.md
cp readme.template README.md

echo "" >> README.md
echo "" >> README.md
echo "\`\`\`bash" >> README.md
./target/debug/modkit --help >> README.md
echo "\`\`\`" >> README.md

subcommands=("adjust-mods" "update-tags" "pileup" "sample-probs" "summary" "motif-bed")
for cmd in "${subcommands[@]}"
do
  echo "" >> README.md
  echo "## ${cmd}" >> README.md
  echo "\`\`\`bash" >> README.md
  ./target/debug/modkit $cmd --help >> README.md
  echo "\`\`\`" >> README.md
done

