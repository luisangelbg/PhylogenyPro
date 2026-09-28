#!/bin/sh
# Rebuilds js/examples.js from the files in data/.
# Run from the project root:  sh data/hacer_examples.sh > js/examples.js
fasta_to_js() {   # $1 = file, $2 = key
  awk -v key="$2" '
    BEGIN { printf "  %s: {\n", key }
    /^>/ { if (n) printf "\",\n"; printf "    \"%s\": \"", substr($0,2); n=1; next }
    { printf "%s", $0 }
    END { if (n) printf "\",\n"; printf "  },\n" }
  ' "$1"
}
echo "/* PhylogenyPro — molecular phylogenetics in the browser."
echo "   Copyright (C) 2026  Luis Ángel Barrera-Guzmán"
echo ""
echo "   This program is free software: you can redistribute it and/or modify it under"
echo "   the terms of the GNU General Public License as published by the Free Software"
echo "   Foundation, either version 3 of the License, or (at your option) any later"
echo "   version. It is distributed in the hope that it will be useful, but WITHOUT ANY"
echo "   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR"
echo "   A PARTICULAR PURPOSE. See the GNU General Public License, in the file LICENSE"
echo "   at the root of this program, or <https://www.gnu.org/licenses/>. */"
echo ""
echo "/* PhylogenyPro — example data sets, embedded so they load without a server."
echo "   Rebuild with:  sh data/hacer_examples.sh > js/examples.js"
echo "   The sequences are SIMULATED on a known tree (data/generar_ejemplos.R); the"
echo "   species names are real, the sequences are not. */"
echo "const EXAMPLES = {"
echo "seqs: {"
fasta_to_js data/ejemplo_bursera/rbcL.fasta rbcL
fasta_to_js data/ejemplo_bursera/trnLF.fasta trnLF
fasta_to_js data/ejemplo_bursera/ITS.fasta ITS
fasta_to_js data/ejemplo_bursera/concatenado.fasta concatenado
fasta_to_js data/ejemplo_proteina.fasta proteina
echo "},"
echo "trees: {"
printf '  bursera: "%s",\n' "$(cat data/ejemplo_bursera/true_tree.nwk)"
printf '  proteina: "%s",\n' "$(cat data/ejemplo_proteina_true_tree.nwk)"
printf '  morfologia: "%s",\n' "$(cat data/ejemplo_morfologia_true_tree.nwk)"
echo "},"
echo "morfologia: \`"
cat data/ejemplo_morfologia.csv
echo "\`,"
echo "particiones: \`"
cat data/ejemplo_bursera/particiones.txt
echo "\`,"
echo "};"
echo "window.EXAMPLES = EXAMPLES;"
