#!/bin/bash

# Test variant decomposition
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Variant decomposition test ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
			 	${DIR}/../../example/decomp/decomp_test.truth.vcf.gz \
			 	${DIR}/../../example/decomp/decomp_test.query.vcf.gz \
			 	-f ${DIR}/../../example/decomp/decomp_test.conf.bed.gz \
			 	-o ${TMP_OUT} \
                --preprocess-truth \
			 	-X -V \
			 	--force-interactive  #  --verbose

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/decomp/expected.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../../example/decomp/expected.summary.csv --schema ga4gh-summary-v2
if [[ $? != 0 ]]; then
	echo "Summary differs! ${TMP_OUT}.summary.csv ${DIR}/../../example/decomp/expected.summary.csv"
	exit 1
fi

gunzip -c ${TMP_OUT}.vcf.gz > ${TMP_OUT}.vcf
diff -I ^# ${TMP_OUT}.vcf ${DIR}/../../example/decomp/expected.vcf
if [[ $? != 0 ]]; then
	echo "Variants differ! diff ${TMP_OUT}.vcf ${DIR}/../../example/decomp/expected.vcf"
	exit 1
fi

rm -rf ${TMP_OUT}.*

echo "Variant decomposition test was successful."
