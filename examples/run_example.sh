



DATA_DIR="$(dirname "$0")/data"
SCRIPT_DIR="$(dirname "$0")"
source "${SCRIPT_DIR}/../bin/activate_cphasing"


echo "Clean up previous test outputs..."
rm -rf cphasing_output_porec cphasing_output_hic

echo "👉 Testing Pore-C Pipeline..."
cphasing pipeline \
    -f "${DATA_DIR}/contigs.fasta.gz" \
    -pcd "${DATA_DIR}/porec_reads.fasta.gz" \
    -t 10 \
    -hcr \
    -p AAGCTT \
    -n 5:4 \
    -o "${SCRIPT_DIR}/cphasing_output_porec"

echo "✅ Pore-C Pipeline test completed successfully!"


echo "👉 Testing Hi-C Pipeline..."
cphasing pipeline \
    -f "${DATA_DIR}/contigs.fasta.gz" \
    -hic1 "${DATA_DIR}/hic_R1.fasta.gz" \
    -hic2 "${DATA_DIR}/hic_R2.fasta.gz" \
    -t 10 \
    -hcr \
    -p AAGCTT \
    -n 5:4 \
    -o "${SCRIPT_DIR}/cphasing_output_hic"

echo "✅ Hi-C Pipeline test completed successfully!"DATA_DIR="$(dirname "$0")/data"
SCRIPT_DIR="$(dirname "$0")"
source "${SCRIPT_DIR}/../bin/activate_cphasing"


echo "Clean up previous test outputs..."
rm -rf cphasing_output_porec cphasing_output_hic

echo "👉 Testing Pore-C Pipeline..."
if cphasing pipeline \
    -f "${DATA_DIR}/contigs.fasta.gz" \
    -pcd "${DATA_DIR}/porec_reads.fasta.gz" \
    -t 10 \
    -hcr \
    -p AAGCTT \
    -n 5:4 \
    -o "${SCRIPT_DIR}/cphasing_output_porec"; then
    echo "✅ Pore-C Pipeline test completed successfully!"
else
    echo "❌ Error: Pore-C Pipeline test failed!" >&2
    exit 1
fi


echo "👉 Testing Hi-C Pipeline..."
if cphasing pipeline \
    -f "${DATA_DIR}/contigs.fasta.gz" \
    -hic1 "${DATA_DIR}/hic_R1.fasta.gz" \
    -hic2 "${DATA_DIR}/hic_R2.fasta.gz" \
    -t 10 \
    -hcr \
    -p AAGCTT \
    -n 5:4 \
    -o "${SCRIPT_DIR}/cphasing_output_hic"; then
    echo "✅ Hi-C Pipeline test completed successfully!"
else
    echo "❌ Error: Hi-C Pipeline test failed!" >&2
    exit 1
fi