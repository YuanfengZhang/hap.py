#!/bin/bash

# Test invalid path traversals
#

set +e

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
. ${DIR}/detect_vars.sh

echo "Test for path traversals: ${HCVERSION} from ${HCDIR}"

TMP_OUT=`mktemp -t happy.XXXXXXXXXX`

# run hap.py
${PYTHON} ${HCDIR}/hap.py \
				${DIR}/../data/pathtraversal/test.vcf \
				${DIR}/../data/pathtraversal/test2.vcf \
				-o ${TMP_OUT} \
				-X --reference ${DIR}/../data/pathtraversal/test.fa -l chrQ \
				--force-interactive

if [[ $? != 0 ]]; then
	echo "hap.py failed!"
	exit 1
fi

mkdir -p "${TMPDIR}/summary-legacy-evidence" || exit 1

${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/pathtraversal/expected.summary.csv > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.log" 2>&1
LEGACY_SUMMARY_STATUS=$?
printf '%s\n' "$LEGACY_SUMMARY_STATUS" > "${TMPDIR}/summary-legacy-evidence/$(basename "${TMP_OUT}.summary.csv").legacy-compare.exit-code"
${PYTHON} ${DIR}/compare_summaries.py ${TMP_OUT}.summary.csv ${DIR}/../data/pathtraversal/expected.summary.csv --schema ga4gh-summary-v2
if [[ $? != 0 ]]; then
	echo "Summary differs! ${TMP_OUT}.summary.csv ${DIR}/../data/pathtraversal/expected.summary.csv"
	exit 1
else
    echo "Path traversal test successful"
    rm -rf ${TMP_OUT}.*
fi

