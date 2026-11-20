#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
quickld="$root/quickLD"
loci=${1:-10}
haplotypes=${2:-10}
seed=${3:-17}
output_dir=${4:-"$root/quickld_generated"}

# qLD-parse-VCF returns 1 after successfully writing its parts. Treat that
# documented-in-practice status as success while preserving real failures.
run_quickld() {
    set +e
    "$@"
    local status=$?
    set -e
    if (( status != 0 && status != 1 )); then
        return "$status"
    fi
}

for executable in "$quickld/bin/qLD-parse-VCF" \
                  "$quickld/bin/qLD-parse-2MDF" \
                  "$quickld/bin/qLD-compute"; do
    if [[ ! -x "$executable" ]]; then
        echo "quickLD is not built. Run: (cd \"$quickld\" && ./install_blis && ./config install blis)" >&2
        exit 2
    fi
done

make -C "$root" generate_hls_vcf
mkdir -p "$output_dir"
vcf="$output_dir/hls_generated.vcf"
parsed="$output_dir/parsed_vcf"
mdf="$output_dir/mdf"
report="$output_dir/quickld_report.txt"

# Every HLS bit is represented as one haploid VCF sample, so this works for
# all supported HLS sample counts, including odd values.
"$root/generate_hls_vcf" "$loci" "$haplotypes" "$seed" "$vcf"
run_quickld "$quickld/bin/qLD-parse-VCF" -input "$vcf" -output "$parsed" -size 1 -chrom HLS
run_quickld "$quickld/bin/qLD-parse-2MDF" -input "$parsed" -output "$mdf" -ploidy haploid -seed "$seed"
run_quickld "$quickld/bin/qLD-compute" -input "$mdf" -output "$report" -ploidy haploid \
    -mdf -blis -threads 1 -r2limit 0

echo "quickLD report: $report"
