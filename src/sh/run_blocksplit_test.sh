#!/bin/bash

##############################################################
# Test setup
##############################################################

set -e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

. ${DIR}/detect_vars.sh

##############################################################
# Test blocksplit
##############################################################


echo "Running blocksplit test"

TF_r=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

VCF1=${DIR}/../../example/happy/PG_NA12878_hg38-chr21.vcf.gz
VCF2=${DIR}/../../example/happy/NA12878-GATK3-chr21.vcf.gz

${HCDIR}/blocksplit $VCF1 $VCF2 \
	-o ${TF_r} -l chr21 -w 10000

${PYTHON} ${DIR}/../python/ovc.py ${TF_r}
if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- overlaps found. You can inspect ${TF_r} for the failed result."
	exit 1
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

${HCDIR}/bcftools view -H $VCF1 chr21 > ${TF_X1}
${HCDIR}/bcftools view -H $VCF1 chr21 -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi

TF_X1=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1
TF_X2=`mktemp -t blocksplit.XXXXXXXXXX` || exit 1

${HCDIR}/bcftools view -H $VCF2 chr21 > ${TF_X1}
${HCDIR}/bcftools view -H $VCF2 chr21 -R ${TF_r} > ${TF_X2}

diff ${TF_X1} ${TF_X2}

if [ $? -ne 0 ]; then
	echo "blocksplit test FAILED -- doesn't cover all variants. You can inspect ${TF_r} / ${TF_X1} / ${TF_X2} for the failed result."
	exit 1
else
	rm -f ${TF_X1} ${TF_X2}
fi

echo "Blocksplit test was successful"

# Region/reader cleanup regression. Never remove the source fixture's index.
REGRESSION_DIR=$(mktemp -d "${TMPDIR}/blocksplit-regions.XXXXXXXX")
trap 'rm -rf "$REGRESSION_DIR"' EXIT
FORMAT_VCF="${DIR}/../data/test_formats.vcf.gz"
"${HCDIR}/blocksplit" "$FORMAT_VCF" -l chr1 -o "$REGRESSION_DIR/present.bed" 2> "$REGRESSION_DIR/present.stderr"
test -s "$REGRESSION_DIR/present.bed"
"${HCDIR}/blocksplit" "$FORMAT_VCF" -l absent_custom_contig -o "$REGRESSION_DIR/absent.bed" 2> "$REGRESSION_DIR/absent.stderr"
test ! -s "$REGRESSION_DIR/absent.bed"
cp "$FORMAT_VCF" "$REGRESSION_DIR/unindexed.vcf.gz"
status=0
"${HCDIR}/blocksplit" "$REGRESSION_DIR/unindexed.vcf.gz" -l chr1 -o "$REGRESSION_DIR/unindexed.bed" 2> "$REGRESSION_DIR/unindexed.stderr" || status=$?
test "$status" -eq 1
if grep -E 'AddressSanitizer|LeakSanitizer|UndefinedBehaviorSanitizer|runtime error:' "$REGRESSION_DIR/"*.stderr; then
    echo "blocksplit region regression emitted a sanitizer error" >&2
    exit 1
fi
echo "Blocksplit region and unindexed-input regressions passed"
